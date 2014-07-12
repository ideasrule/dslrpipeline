#!/usr/bin/python -u
""" Implements the masterphotref action, which generates a master photometric
reference. """

from MagnitudeFitting import get_config, \
                                              SinglePhotRef, \
                                              start_grcollect, \
                                              all_linear_terms,\
                                              read_master_catalogue,\
                                              start_magfit, log_magfit,\
                                              mag_from_mag_per_minute,\
                                              required_sphotref_keys,\
                                              patch_header_keywords,\
                                              read_config_file,\
                                              output_done,\
                                              get_phot_type
from HATpipepy.Common.DBUtil import raw_db_table, get_db, column_in_list
from HATpipepy.HATNet import Stations as HATNetStations
from HATpipepy.HATSouth import Stations as HATSouthStations
from HATpipepy.Common import HATUtil
from HATpipepy.Common.HATUtil import get_night_str, execute_pipe, \
                                     parse_hat_id, parse_frame_name,\
                                     parse_key_string
from HATpipepy.Common.FitsUtil import get_keywords, get_image_vars
from HATpipepy.Common.CatUtil import Cat2MASS, is_AAA_quality_2MASS, \
                                     JmK_color_cut_2MASS
from HATpipepy.Common.AstromUtil import trans_file_center, FilterSourceList
from HATpipepy.Common.UglyDB import db_median
from HATpipepy.Common.CalibStatus import object_status
from HATpipepy.Common.CalibRecover import RecoveryInfo
from HATpipepy.Common.FixCrossmountSymlinks import fix_symlinks_for
from HATpipepy.Common import Error
from optparse import OptionParser
from os.path import basename, splitext, exists, dirname, expanduser
from os import makedirs, rename, chmod, kill
from os.path import join as join_paths
from math import ceil
from subprocess import Popen, PIPE, call
from scipy import median, nan, isnan, array, zeros, dot, shape, sqrt, \
                  delete, log10
from scipy.optimize import leastsq
import logging
import sys
from signal import SIGTERM

SVNinfo=dict(Revision=eval('$Revision: 2617 $'[10:-1]),
             Date='$Date: 2013-08-12 17:52:49 -0400 (Mon, 12 Aug 2013) $'[7:-1],
             Author='$Author: kpenev $'[9:-1],
             URL='$URL: https://hat.astro.princeton.edu/svn/HATreduc/HATpipe/source/HATpipepy/Actions/do_masterphotref.py $'[6:-1])
log_extra=dict(filever=SVNinfo['Revision'])

def get_work(db, required_hdr_keys, restrict_to_fields,
             logger) :
    """ Figures out the frames for which single magnitude fitting is pending.
    
    Returns:
        one nested dictionary per frame type (object, focus, cal) indexed by
        project_id, field, camera position and single photometric reference
        id, containing a list of tuples (single photometric reference header,
        a list of header dictionaries for frames which should be combined to
        generate a master photometry reference). """

    def get_frame_list(frame_type) :
        """ Returns a list of candidate frames to do single photometric
        refence fitting for. """

        frame_tbl=raw_db_table(frame_type)
        frame_list=db("SELECT `F`.`project_id`, `F`.`object`, `F`.`cmpos`, "
                      "`F`.`station_id`, `F`.`fnum`, `F`.`night` FROM `"+
                      frame_tbl+"` AS `F` JOIN `SinglePhotRef` AS `S` ON "
                      "(`F`.`project_id`=`S`.`project_id` AND "
                      "`F`.`object`=`S`.`object` AND "
                      "`F`.`cmpos`=`S`.`cmpos` AND `S`.`imagetype`=%s) "
                      "LEFT JOIN `MasterPhotRef` AS"
                      " `M` ON (`S`.`project_id`=`M`.`project_id` AND "
                      "`S`.`sphotref_id`=`M`.`sphotref_id`) WHERE"
                      " `S`.`enabled`='1' AND `F`.`calib_status`>=%s AND "
                      "(`M`.`copies` IS NULL OR `M`.`copies`='') GROUP BY "
                      "`F`.`project_id`, `F`.`object`, `F`.`cmpos`, "
                      "`F`.`station_id`, `F`.`fnum`", 
                      (frame_type, object_status['raw_photometry']),
                      no_simplify=True)
        logger.debug('database work querry returned %d frames.'%(
            len(frame_list) if frame_list is not None else 0),
            extra=log_extra)
        if frame_list is None : return []
        else : return frame_list

    def process(frame_record) :
        """ Returns True if the given frame should be processed. """

        if not restrict_to_fields : return True
        for allowed_field, allowed_cams in restrict_to_fields :
            if (frame_record[1]==allowed_field and 
                (allowed_cams==() or frame_record[2] in allowed_cams)) :
                return True
        return False

    result=[]
    frame_keys=dict()
    for key_str in db('SELECT `keys` FROM `SinglePhotRef`', 
                      no_simplify=True) :
        frame_keys.update(parse_key_string(key_str[0]))
    for k in ['PROJID', 'OBJECT', 'CMPOS', 'STID', 'FNUM', 'NIGHT'] :
        if k in frame_keys : del frame_keys[k]
    existing_masters=db("SELECT `project_id`, `sphotref_id` FROM "
                        "`MasterPhotRef` WHERE `copies` IS NOT NULL AND "
                        "`copies`!=''", (), no_simplify=True)
    if existing_masters is None : existing_masters=[]
    for frame_type in ['object', 'focus', 'cal'] :
        logger.debug('Getting work for %s frames.'%frame_type,
                     extra=log_extra)
        frame_list=filter(process, get_frame_list(frame_type))
        work=dict()
        current_project=None
        current_field=None
        current_cmpos=None
        dest=None
        for f in frame_list :
            project=f[0]
            field=f[1]
            cmpos=f[2]
            if (project!=current_project or field!=current_field or
                cmpos!=current_cmpos) :
                get_single_photref=SinglePhotRef(db, project, field, cmpos,
                                                 frame_type,
                                                 required_hdr_keys)
                if project!=current_project :
                    work[project]=dict()
                current_project=project
                current_field=field
                current_cmpos=cmpos
            frame_hdr=dict(PROJID=project, OBJECT=field, CMPOS=cmpos,
                           STID=f[3], FNUM=f[4], NIGHT=get_night_str(f[5]),
                           IMAGETYP=frame_type)
            frame_hdr.update(get_keywords(db.red_frame(
                f[3], frame_hdr['NIGHT'], f[4], cmpos), frame_keys))
            spr_header=get_single_photref(frame_hdr)
            if ((spr_header is None) or 
                ((spr_header['PROJID'], spr_header['SPRID']) in
                 existing_masters)) : continue
            sphotref_id=spr_header['SPRID']
            if sphotref_id not in work[project] :
                work[project][sphotref_id]=(spr_header, [])
            work[project][sphotref_id][1].append(frame_hdr)
        result.append(work)
    return result

class dummy_class() : pass

class TransFname :
    """ Returns the transformation file corresponding to the given header.
    """

    def __init__(self, db) :
        """ Reads the necessary templates from the database. """

        self.__db=db
        self.__astrom_template=db.get_name_template(db.station, 'astrometry',
                                                    'outpath')[0]
    
    def __call__(self, header) :
        """ Returns the transformation file corresponding to the given
        header. """

        fits, packed=self.__db.red_frame(header['STID'], header['NIGHT'],
                                         header['FNUM'], header['CMPOS'], 
                                         add_packed=True)
        if packed : fits=splitext(fits)[0]
        trans=(splitext(basename(fits))[0]+
               self.__db.extensions['astrom']['trans'])
        return join_paths(self.__astrom_template.substitute(header), trans)

def catalogue_fov(db, frames, xres, yres, config, trans_fname, logger,
                  recover) :
    """ Returns a center (ra,dec) and a field of view (degrees) which cover
    all frames correposponding to the given transformation files. """

    def get_querry_center(trans_files, max_outside) :
        """ Returns the optimal RA and DEC around which the catalogue querry
        should be done. """

        join=''
        obj_status=''
        args=[]
        for i, imagetype in enumerate(['object', 'focus', 'cal']) :
            short='`F%d`'%i
            join+=(' LEFT JOIN `'+raw_db_table(imagetype)+'` AS '+short+
                   ' ON (`P`.`station_id`='+short+'.`station_id` AND '
                   '`P`.`fnum`='+short+'.`fnum` AND `P`.`cmpos`='+short+
                   '.`cmpos`)')
        limit=(max_outside/2, len(trans_files)-max_outside/2)
        querry=('FROM `AstromPointing` AS `P` '+join+' WHERE '
                '`P`.`station_id`=%s AND `P`.`fnum`=%s AND `P`.`cmpos`=%s ')
        ra_dec=map(lambda arg: db('SELECT 15.0*`ra`, `dec` '+querry, arg),
                   map(lambda h: (h['STID'], h['FNUM'], h['CMPOS']), frames))
        ra=sorted(map(lambda rd: rd[0], ra_dec))
        dec=sorted(map(lambda rd: rd[1], ra_dec))
        assert len(ra)==len(trans_files)
        assert len(ra)==len(dec)
        ra=0.5*(ra[limit[0]]+ra[limit[1]])
        dec=0.5*(dec[limit[0]]+dec[limit[1]])
        return (round(ra, config.mcat.round_pointing),
                round(dec, config.mcat.round_pointing))

    def get_querry_size(trans_files, ra_global, dec_global, max_outside) :
        """ Returns the field of view that should be fetched from the
        catalogue in order to cover all the frames corresponding to the given
        transformation files. """

        work=map(lambda f : (f,)+trans_file_center(f), trans_files)
        var_arg_fmt='arc,ra=%s,dec=%s,degrees'
        trans_pipe=[['grtrans', '--input', '-', '--col-xy', '1,2',
                     '--input-transformation', None, '--reverse'],
                    ['grtrans', '--col-pixel', '1,2', '--wcs', None],
                    ['grtrans', '--col-radec', '1,2', '--wcs', 
                     var_arg_fmt%(repr(ra_global), repr(dec_global))]]
        pipe_input="0 0\n0 %d\n%d 0\n%d %d\n"%(yres, xres, xres, yres)
        fov=0
        current_outliers=set()
        current_min_outlier=0
        for hdr, (fname, ra, dec) in zip(frames, work) :
            trans_pipe[1][-1]=var_arg_fmt%(repr(ra), repr(dec))
            trans_pipe[0][6]=fname
            running=[Popen(trans_pipe[0], stdin=PIPE, stdout=PIPE)]
            for i in range(1, len(trans_pipe)) :
                running.append(Popen(trans_pipe[i],stdin=running[i-1].stdout,
                                     stdout=PIPE))
            running[0].stdin.write(pipe_input)
            old_fov=fov
            for i in range(4) :
                l=running[-1].stdout.readline().strip()
                xi, eta=map(lambda v: abs(eval(v)), l.split())
                dist=max(xi, eta)
                if (max_outside!=0 and dist>current_min_outlier) :
                    current_outliers.add(dist)
                    if len(current_outliers)>max_outside :
                        current_outliers.remove(current_min_outlier)
                        fov=max(current_min_outlier, fov)
                        current_min_outlier=min(current_outliers)
                    elif current_min_outlier<dist :
                        current_min_outlier=dist
                else : fov=max(dist, fov)
            for r in running : 
                kill(r.pid, SIGTERM) 
                r.wait()
            if 2.0*fov>config.mcat.fov_alarm :
                fov=old_fov
                logger.warning("The transformation file '%s' places at "
                               "least one of the corners of the "
                               "corresponding frame further than %f "
                               "degrees away from the catalogue querry "
                               "center. Changing calibration status to "
                               "badpointing."%
                               (fname, config.mcat.fov_alarm),
                               extra=log_extra)
                db_key=(hdr['STID'], hdr['FNUM'], hdr['CMPOS'])
                for imtype in ['object', 'focus', 'cal'] :
                    db('UPDATE `'+raw_db_table(imtype)+'` SET '
                       '`calib_status`=%s WHERE `station_id`=%s AND '
                       '`fnum`=%s AND `cmpos`=%s',
                       (object_status['bad_pointing'],)+db_key)
        return 2.0*fov

    trans_files=map(trans_fname, frames)
    max_outside=int(config.mphotref.min_measurements-1 if 
                    config.mphotref.min_measurements>=1
                    else (config.mphotref.min_measurements*len(trans_files)))
    ra, dec=get_querry_center(trans_files, max_outside)
    fov=ceil(10.0*config.mcat.fov_safety_fac*get_querry_size(
        trans_files, ra, dec, max_outside))/10.0
    return ra, dec, fov

def catalogue_exists(db, project, sphotref_id) :
    """ Returns a pair of boolean values indicating whether the full and
    sparse, respectively, catalogues exist for the given field and
    camera position. """

    return bool(db('SELECT * FROM `MasterCatalogue` WHERE `project_id`=%s '
                   'AND `sphotref_id`=%s AND `expired`=%s',
                   (project, sphotref_id, 0)))

def generate_master_catalogue(db, config, sph_ref_hdr,
                              fov, recover, logger, insert_in_database=False,
                              mphotref_ver=None, rawphot_ver=None) :
    """ Returns a master catalogue satisfying the given parameters - a
    dictionary indexed by (field, source) identifiers contaning dictionaries
    of relevant 2MASS information. """

    def update_db() :
        """ Updates the database to indicate that a new master catalogue
        querry was generated. """

        if insert_in_database :
            db('INSERT INTO `MasterCatalogue` SET `project_id`=%s, '
               '`sphotref_id`=%s, `object`=%s, `cmpos`=%s, '
               '`mphotref_ver`=%s, `rawphot_ver`=%s, `ra`=%s, `dec`=%s, '
               '`fieldofview`=%s, `copies`=%s',
               (sph_ref_hdr['PROJID'], sph_ref_hdr['SPRID'], 
                sph_ref_hdr['OBJECT'], sph_ref_hdr['CMPOS'], mphotref_ver,
                rawphot_ver, fov[0], fov[1], fov[2], str(db.station)))

    master_fname=config.mcat.template.substitute(sph_ref_hdr)
    master_dir=dirname(master_fname)
    if not exists(master_dir) : makedirs(master_dir, 0775)
    cat_querry=dict(output=config.mcat.columns,
                    ra=fov[0], dec=fov[1], fov=fov[2],
                    xi_eta_reference=(sph_ref_hdr['NRAC'], 
                                      sph_ref_hdr['NDECC']),
                    filter=sph_ref_hdr['FILTERS'],
                    bright=config.mcat.bright_mag,
                    faint=config.mcat.faint_mag)
    recover.start_section()
    proxy_fname=Cat2MASS(db, recover, logger)(**cat_querry)
    recover.drop_section()
    if call(['mv', proxy_fname, master_fname]) :
        raise Error.External('Moving %s -> %s failed.'%(proxy_fname,
                                                        master_fname))
    chmod(master_fname, 0775)
    update_db()

def master_catalogue(db, sphotref_header, frames, config, recover, logger) :
    """ Generates the master catalogue files necessary to do the magnitude
    fitting indicated in work. """

    trans_fname=TransFname(db)
    project=sphotref_header['PROJID']
    sphotref_id=sphotref_header['SPRID']
    sph_trans=trans_fname(sphotref_header)
    catalogue_fname=config.mcat.template.substitute(sphotref_header)

    if catalogue_exists(db, project, sphotref_id) :
        logger.info("Catalogue '%s' already exists for project %d, single "
                    "photometric reference %d, not re-generating."%
                    (catalogue_fname, project, sphotref_id), extra=log_extra)
    else :
        logger.info("Generating Catalogue '%s' for project %d, single "
                    "photometric reference %d."%
                    (catalogue_fname, project, sphotref_id), extra=log_extra)
        cat_fov=catalogue_fov(db, frames, sphotref_header['NAXIS1'],
                              sphotref_header['NAXIS2'], 
                              config, trans_fname, logger, None)
        logger.info('Catalogue field of view for project %d, single '
                    'photometric reference %d: %f'%
                    (project, sphotref_id, cat_fov[2]), extra=log_extra)
        generate_master_catalogue(db, config, sphotref_header,
                                  cat_fov, recover, logger, True,
                                  config.mphotref.version,
                                  config.rawphot_ver)
        logger.info('Generated master catalogue for field %s, '
                    'camera %d'%(sphotref_header['OBJECT'],
                                 sphotref_header['CMPOS']), extra=log_extra)

def process_magfit_stat(stat_columns, sph_header, config, db,
                        num_input_files, magfit_ver) :
    """ Generates a master photometric reference from the statistics file of
    the single photometric reference magnitude fitting. """

    def get_med_count(stat_fname) :
        """ Returns the median number of observations for the sources in the
        given statistics file. """

        if 'rcount' in stat_columns :
            count_col=stat_columns.index('rcount')
        else : count_col=stat_columns.index('count')
        stat_file=open(stat_fname, 'r')
        med=median(map(lambda l: eval(l.split()[count_col]), stat_file))
        stat_file.close()
        return med

    def get_stat(stat_fname, master_cat, min_counts) :
        """ Returns a statistics structure consisting of a list of
        dictionaries (one for each aperture) with keys source - HAT-.., 
        M - list of average instrumental magnitudes, scatter - list
        of scatter estimates of the given magnitudes. """
        stat_file=open(stat_fname, 'r')            
        mag_letter=sph_header['FILTERS']
        stat=[{'source':[], 'full count':[], 'rej count':[], mag_letter:[], 
               'scatter':[], 'xi':[], 'eta':[], 'scatter2':[]} 
              for i in range(config.num_apertures)]
        full_count_col=stat_columns.index('count')
        id_col=stat_columns.index('id')
        if 'rcount' in stat_columns :
            rej_count_col=stat_columns.index('rcount')
            med_col=stat_columns.index('rmedian')
            scatter_col=stat_columns.index('rmedianmeddev')
            scatter2_col=stat_columns.index('rmediandev')
        else : 
            rej_count_col=full_count_col
            med_col=stat_columns.index('median')
            scatter_col=stat_columns.index('medianmeddev')
            scatter2_col=stat_columns.index('mediandev')
        col_per_ap=len(stat_columns)-1
        saturation_mag=mag_from_mag_per_minute(
            config.mphotref.rms_fit_bright_mag_min, sph_header['EXPTIME'])
        for l in stat_file :
            entries=l.split()
            for ap in range(config.num_apertures) :
                rej_count=eval(entries[rej_count_col+ap*col_per_ap])
                if rej_count<min_counts : continue
                hatid=entries[id_col]
                if (hatid not in master_cat or
                    master_cat[hatid][mag_letter]<saturation_mag) : continue
                mag=eval(entries[med_col+ap*col_per_ap])
                scatter=eval(entries[scatter_col+ap*col_per_ap])
                if isnan(mag) or isnan(scatter) or scatter<=0 : continue
                stat[ap]['source'].append(hatid)
                stat[ap]['full count'].append(eval(entries[full_count_col]))
                stat[ap]['rej count'].append(rej_count)
                stat[ap][mag_letter].append(mag)
                stat[ap]['scatter'].append(scatter)
                stat[ap]['scatter2'].append(eval(entries[scatter2_col]))
                for col in ['xi', 'eta'] :
                    stat[ap][col].append(master_cat[hatid][col])
        return stat

    def rejected_indices(residuals) :
        """ Returns a list of indices for which the residuals are large
        enough to be considered outliers. """

        fit_res=sqrt(config.mphotref.rms_fit_err_avg(pow(residuals, 2)))
        max_res=fit_res*config.mphotref.rms_fit_rej_lvl
        return filter(lambda i: abs(residuals[i])>max_res, 
                      range(len(residuals)-1,-1,-1)), fit_res

    def fit_to_db(coefficients, ap_ind, fit_res, start_src_count,
                  final_src_count) :

        db_columns=all_linear_terms(config.mphotref.rms_fit_param)[0]
        db('REPLACE INTO `MasterPhotRef` (`project_id`, `sphotref_id`, '
           '`object`, `cmpos`, `station_id`, `aperture`, `imagetype`, '
           '`copies`, `magfitver`, `mphotrefcfgver`, `input_src`, '
           '`non_rej_src`, `fit_residual`, `'+
           '`, `'.join(db_columns)+'`) VALUES (%s'+
           ', %s'*(len(db_columns)+12)+')', 
           (sph_header['PROJID'], sph_header['SPRID'], sph_header['OBJECT'],
            sph_header['CMPOS'], sph_header['STID'], ap_ind,
            sph_header['IMAGETYP'], '', magfit_ver, config.mphotref.version,
            start_src_count, final_src_count, fit_res)+tuple(coefficients))

    def fit_scatter(ap_ind, stat) :
        """ Derives a fit of the scatter from the given singe aperture
        statistics dictionary. """

        result=[]
        all_deriv=all_linear_terms(config.mphotref.rms_fit_param, stat)[0]
        derivatives=all_deriv[:]
        num_free_coef=len(derivatives)
        alllogscatter=log10(array(stat['scatter']))
        logscatter=alllogscatter
        start_source_count=len(logscatter)
        error_func=lambda coef: dot(coef, derivatives)-logscatter
        deriv_func=lambda coef: derivatives
        initial_guess=zeros(num_free_coef)
        for rej_iter in range(config.mphotref.rej_iterations) :
            if len(logscatter)<num_free_coef : 
                print sph_header
                print len(logscatter)
                raise Error.Numeric("Rejecting outliers resulted in less "
                                    "points than fit parameters when "
                                    "fitting log10(magnitude scatter) "
                                    "for aperture %d of field %s, camera %s,"
                                    " %s frames."%
                                    (ap_ind, sph_header['OBJECT'],
                                     sph_header['CMPOS'],
                                     sph_header['IMAGETYP']))
            coefficients, covariance, info_dict, msg, status=leastsq(
                error_func, initial_guess, Dfun=deriv_func, col_deriv=1,
                full_output=1)
            if status not in [1, 2, 3, 4] :
                raise Error.Numeric("Linear least squares fitting of the "
                                    "scatter of the single photometry "
                                    "magnitudes for aperture %d of field %s,"
                                    " camera %s, %s frames failed: %s"%
                                    (ap_ind, sph_header['OBJECT'],
                                     sph_header['CMPOS'],
                                     sph_header['IMAGETYP'], msg))
            bad_ind, fit_res=rejected_indices(info_dict['fvec'])
            if not bad_ind : break
            if rej_iter==config.mphotref.rej_iterations-1 : break
            derivatives=map(lambda d : delete(d, bad_ind), derivatives)
            logscatter=delete(logscatter, bad_ind)
        if db is not None :
            fit_to_db(coefficients, ap_ind, fit_res, start_source_count,
                      len(logscatter))
        return alllogscatter-dot(coefficients, all_deriv), fit_res

    def generate_masters(stat) :
        """ Generates the final master file from the statistics structure.
        """

        fmt=' '.join(['%s', '%6d', '%6d', '%12.4f', '%12.5f', '%12.7f', 
                      '%12.5f'])+'\n'
        colname_fmt='#'+' '.join(['%14s', '%6s', '%6s']+['%12s']*4)+'\n'
        mag_letter=sph_header['FILTERS']
        for ap_ind, ap_stat in enumerate(stat) :
            master=open(config.mphotref.template.substitute(
                APERTURE=ap_ind, **sph_header), 'w')
            master.write(colname_fmt%(('ID', 'full_count', 'rej_count', 
                                       'mag', 'scatter', 'scatter excess',
                                       'scatter2')))
            master.write(colname_fmt%tuple('['+str(i+1)+']' 
                                           for i in range(7)))
            excess=ap_stat['scatter excess'][:]
            excess.sort()
            max_scatter_excess=min(
                config.mphotref.max_rms_above_fit*ap_stat['fit residual'],
                excess[-int(config.mphotref.max_rms_quantile*len(excess))])
            records=zip(ap_stat['source'], 
                        ap_stat['full count'], ap_stat['rej count'],
                        ap_stat[mag_letter], ap_stat['scatter'],
                        ap_stat['scatter excess'], ap_stat['scatter2'])
            records=filter(lambda r: r[-2]<max_scatter_excess, records)
            records.sort()
            map(lambda r: master.write(fmt%r), records)
            master.close()
            if db is not None :
                db("UPDATE `MasterPhotRef` SET `copies`=%s WHERE "
                   "`project_id`=%s AND `sphotref_id`=%s AND `aperture`=%s",
                   (str(db.station), sph_header['PROJID'],
                    sph_header['SPRID'], ap_ind))

    stat_fname=config.magfit.single.stat_template.substitute(sph_header)
    min_counts=config.mphotref.min_meas_med*get_med_count(stat_fname)
    if config.mphotref.min_measurements>=1 : 
        min_counts=max(min_counts, config.mphotref.min_measurements-1)
    else : min_counts=max(min_counts, 
                          config.mphotref.min_measurements*num_input_files)
    master_cat=read_master_catalogue(
        config.mcat.template.substitute(sph_header),
        config.mcat.columns, True)
    stat=get_stat(stat_fname, master_cat, min_counts)
#    print stat
    for i, s in enumerate(stat) : 
        scatter_excess, s['fit residual']=fit_scatter(i, s)
        s['scatter excess']=list(scatter_excess)
    generate_masters(stat)

def do_masterphotref(db, options, logger, recover) :
    """ Generates a master photometry reference for all fields for which a
    single photometric reference has been chosen. """

    required_hdr_keys=required_sphotref_keys(db)
    work=get_work(db, required_hdr_keys, options.restrict_to_fields, logger)
    trans_fname=TransFname(db)
    for imagetype, type_work in zip(['object', 'focus', 'cal'], work) :
        for project, project_work in type_work.iteritems() :
            for sphotref_id, sphotref_work in project_work.iteritems() :
                sph_header, magfit_work=sphotref_work
                config=get_config(sph_header, db)
                master_catalogue(db, sphotref_work[0], sphotref_work[1],
                                 config, recover, logger)
                logger.info("Generating master photometric reference for the"
                            " %s frames of field %s, camera %d from %d "
                            "frames."%(imagetype, sph_header['OBJECT'],
                                       sph_header['CMPOS'],
                                       len(magfit_work)), extra=log_extra)
                sph_trans=trans_fname(sph_header)
                recover.start_section()
                magfit_ver=1
                stat_columns=['id', 'count', 'rcount', 'rmedian', 
                              'rmediandev', 'rmedianmeddev']
                magfit, magfit_ver=start_magfit(project, sphotref_id,
                                                sph_header['OBJECT'],
                                                sph_header['CMPOS'],
                                                imagetype, options, db,
                                                'single')
                tee=Popen(['tee', '/data/scratch/w/G568_grcollect_input'],
                          stdin=magfit.stdout, stdout=PIPE)
                grcollect, stat_columns=start_grcollect(
                    'single', config, sph_header, tee.stdout)
                magfit.stdout.close()
                tee.stdout.close()
                magfit_err=magfit.stderr.read()
                magfit.stderr.close()
                grcollect_out, grcollect_err=grcollect.communicate()
                log_magfit(db, imagetype, sph_header['OBJECT'],
                           sph_header['CMPOS'], logger, 'single')
                if magfit.wait() :
                    raise Error.External('Single photometric reference '
                                         'magnitude fitting failed for the '
                                         '%s frames of field %s, camera %d.'%
                                         (imagetype, sph_header['OBJECT'],
                                          sph_header['CMPOS'])+magfit_err)
                if grcollect.returncode :
                    raise Error.External(
                        'The grcollect command for generating the statistics'
                        ' from the single photometric reference magnitude '
                        'fit failed with return code of %d and error '
                        'message:'%grcollect.returncode+grcollect_out+
                        grcollect_err)
                process_magfit_stat(stat_columns, sph_header, 
                                    config, db, len(magfit_work),
                                    magfit_ver)
                logger.info("Done with the master photometric reference for "
                            "the %s frames of field %s, camera %d from %d "
                            "frames."%(imagetype, sph_header['OBJECT'], 
                                       sph_header['CMPOS'],
                                       len(magfit_work)), extra=log_extra)

def parse_command_line() :
    """ Parses the command line.

    Returns:
        an options structure as produced by OptionParser, but with an
        additional member: phot_list - the list of single reference magnitude
        fitted frames from which a master reference is to be generated.

    Raises:
        Error.CommandLine if parsing the command line fails. """

    parser=OptionParser(usage='%prog <NETWORK> <single photref fits> '
                        '<phot list file>',
                        description='Generates a master photometric '
                        'reference from a list of single reference fitted '
                        'photometry files. Does not touch the database or '
                        'modify the input files in any way.')
    parser.add_option('--log-config', action='store', type='string',
                      dest='log_conf', default='~/CALIBLOG/logging.conf', 
                      help="A file containing the configuration of the "
                      "logger (%default by default). See python logging "
                      "module documentation for the format. The name of the "
                      "logger should be 'LogCalib'.")
    parser.add_option('--config-file', action='store', default='magfit.cfg',
                      help='The name of a file to read various configuration'
                      ' options from. Only used if --manual-frame-list '
                      'options is passed. Default: %default.')
    options, positional_arguments=parser.parse_args()
    if(len(positional_arguments)!=3) :
        parser.print_help()
        raise Error.CommandLine('Exactly three positional arguments '
                                'expected, got %d'%len(positional_arguments))
    network, options.sphotref_fits, options.phot_list=positional_arguments
    if network.upper()=='HATNET' : HATUtil.Stations=HATNetStations
    elif network.upper()=='HATSOUTH' : HATUtil.Stations=HATSouthStations
    else : raise Error.CommandLine("Unrecognized telescope network: "+
                                   network)
    return options

if __name__=='__main__' :
    options=parse_command_line()
    logging.config.fileConfig(expanduser(options.log_conf))
    logger=logging.getLogger("LogCalib")
    logger.info('Outputting %d pre-fitted frames.\n'%len(options.phot_list))
    config=read_config_file(options.config_file)
    required_hdr_keys=required_sphotref_keys()
    sph_header=get_image_vars(options.sphotref_fits, required_hdr_keys)
    phot_list_file=open(options.phot_list, 'r')
    phot_list=[l.strip()  for l in phot_list_file]
    phot_list_file.close()
    patch_header_keywords(sph_header)
    sph_header['PHOT_TYPE']=get_phot_type(phot_list)
    grcollect, stat_columns=start_grcollect('single', config, sph_header,
                                            PIPE)
    output_done(phot_list, 'single', None, config.num_apertures, logger,
                grcollect.stdin)
    grcollect_out, grcollect_err=grcollect.communicate()
    if grcollect.returncode :
        raise Error.External(
            'The grcollect command for generating the statistics'
            ' from the single photometric reference magnitude '
            'fit failed with return code of %d and error '
            'message:'%grcollect.returncode+grcollect_out+
            grcollect_err)
    process_magfit_stat(stat_columns, sph_header, 
                        config, None, len(phot_list), 0)
    logger.info("Done with the master photometric reference for "
                "the %s frames of field %s, camera %d from %d "
                "frames."%(sph_header['IMAGETYP'], sph_header['OBJECT'], 
                           sph_header['CMPOS'],
                           len(phot_list)), extra=log_extra)
    exit(0)



    db=get_db(expanduser('~/CALIBLOG/.calib_db'), logger)
    config=get_config(db)
    work=get_work(db, config.required_hdr_keys)
    trans_fname=TransFname(db)
    field, cmpos=sys.argv[2], eval(sys.argv[3])
    imagetype='object'
    trans_list=map(trans_fname, work[0][field][cmpos][-1])
    print 'deriving field of view for %s camera %d from %d trans files'%(
        field, cmpos, len(trans_list))
#    ra, dec, fov=catalogue_fov(db, trans_list, 4096, 4096, field, cmpos,
#                               config, logger, None)
#    print 'catalogue field: ra=%f, dec=%f, field of view=%f'%(ra,dec,fov)
    ra, dec, fov=None, None, None
    if sys.argv[1]=='catfov' : exit(0)
    sph_header, sph_ver, magfit_work=work[0][field][cmpos]
    sph_trans=trans_fname(sph_header)
    recover=RecoveryInfo('test_recover', logger).new_thread()
#    ra, dec, fov=170.500000, -24.3, 8.1
    generate_master_catalogue(db, field, cmpos, config, sph_trans,
                              sph_header, (ra, dec, fov), 
                              recover, logger)
    if sys.argv[1]=='cat' : exit(0)
    options=dummy_class()
    options.localdb=expanduser('~/CALIBLOG/.calib_db')
    options.max_threads=8
    magfit, magfit_ver=start_magfit(field, cmpos, imagetype, options, db)
    grcollect=start_grcollect('single', config, sph_header,
                              magfit.stdout)
    magfit.stdout.close()
    grcollect_out, grcollect_err=grcollect.communicate()
    log_magfit(imagetype, field, cmpos, logger)
    if magfit.wait() :
        print 'Single photometric reference magnitude fitting failed!'
    if grcollect.returncode :
        print ('The grcollect command for generating the statistics from'
               'the single photometric reference magnitude fit failed '
               'with return code of %d and error message:'%
               grcollect.returncode+grcollect_out+grcollect_err)
#    magfit_ver=1
    process_magfit_stat(sph_header, sph_ver, config, db, len(trans_list),
                        magfit_ver)
