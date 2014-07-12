#!/usr/bin/python2.6 -u

import sys
if sys.version_info[0]>2 or sys.version_info[1]>5 :
    from multiprocessing import Pool, Manager
from scipy import array, ones, dot, zeros, mean, sqrt, isnan, shape, empty,\
                  delete, median, nan, inf, log10
from scipy.optimize import leastsq
from HATpipepy.Common.CalibStatus import object_status
from HATpipepy.Common import HATUtil
from HATpipepy.HATNet import Stations as HATNetStations
from HATpipepy.HATSouth import Stations as HATSouthStations
from HATpipepy.Common.HATUtil import get_night_str, parse_hat_id,\
                                     hat_id_fmt, parse_key_string
from HATpipepy.Common.PhotUtil import raw_phot_cfg_version
from HATpipepy.Common.FitsUtil import get_keywords, get_image_vars
#from HATpipepy.Common.BinPhot import read_fiphot, add_columns
from read_phot import read_fiphot, add_columns
from HATpipepy.Common.DBUtil import get_db, raw_db_table, column_in_list
#from HATpipepy.Common.VerifyOutputFile import check_photometry
from HATpipepy.Common import Error
from traceback import print_exc, format_exception
from optparse import OptionParser
from subprocess import Popen, PIPE
from HATpipepy.Common import CalibLogHandlers
import logging, logging.handlers, logging.config
from HATpipepy.Common import CalibLogHandlers
from string import Template
import os

def check_photometry(*pos_args, **kw_args) : pass

SVNRevision=eval("$Revision: 2277 $"[10:-1])

class dummy_class : pass

def remove_source_ind(phot, source_ind) :
    """ Removes a source from the given photometry file (format should be
    like the output of BinPhot.read_fiphot. """

    for col in phot.keys() :
        if col=='per aperture' :
            for ap in phot[col] :
                for subcol in ap.keys() : del ap[subcol][source_ind]
        elif col!='serial' : del phot[col][source_ind]

class skip_iter :
    """ Iterates over a sequence skipping the sorted sequence of indices
    given at construction. """

    def __init__(self, seq, skip) :
        """ Create the iterator. """

        self.__seq=seq
        self.__skip=iter(skip)
        try : self.__next_skip=self.__skip.next()
        except StopIteration : self.__next_skip=-1
        self.__ind=0

    def __iter__(self) : return self

    def next(self) :
        while self.__ind==self.__next_skip :
            self.__ind+=1
            try : self.__next_skip=self.__skip.next()
            except StopIteration : self.__next_skip=-1
        if self.__ind>=len(self.__seq) : raise StopIteration()
        res=self.__seq[self.__ind]
        self.__ind+=1
        return res

def skip_array(seq, skip) :
    """ Converts the given 1D sequence to a 1D scipy array skipping the given
    indices. """

    if not skip : return array(seq)
    skip_iter=iter(skip)
    try : next_skip=skip_iter.next()
    except StopIteration : next_skip=-1
    res_size=len(seq)-len(skip)
    res=empty(res_size)
    seq_ind=0
    for i in range(res_size) :
        while seq_ind==next_skip : 
            seq_ind+=1
            try : next_skip=skip_iter.next()
            except StopIteration : next_skip=-1
        res.itemset(i,seq[seq_ind])
        seq_ind+=1
    return res

def no_fit_indices(phot, max_mag_err) :
    """ Returns a sorted list of all indices for which at least one of the
    columns in phot is nan. """

    ind_list=range(len(phot['source']))
    return [filter(
        lambda i: ((ap['mag err'][i]>max_mag_err) or isnan(ap['mag'][i]) or
                   isnan(ap['mag err'][i]) or isnan(ap['ref mag'][i]) or
                   isnan(ap['ref mag err'][i])), ind_list) 
        for ap in phot['per aperture']]

def linear_terms_subset(var, xi, eta, max_var_order, max_spatial_order,
                        term_type) :
    """ Returns a list of 
        1 numpy arrays containing terms of independent variables that
          participate in the fit. (term_type='arr')
        2 strings containing the names that should be used for the
          coefficients, either as database column names or as arguments
          to lfit (term_type='coef')
        3 strings containing expressions representing the corresponding
          term in a format that lfit understands. (term_type='expr')
    for terms of a polynomial involving the given variable with up to
    var_order order terms in the variable and up to spatial_order order
    terms in xi and eta.

    The value of var should be a list of 
        * numpy arrays in case 1 
        * variable names in cases 2 and 3. """
    
    def next_power_set(power_set) :
        """ Given a sequence of orders identifying a term in a polynomial
        it modifies it to give the next term of the same total order,
        returning True if there are yet more terms and False if not. """

        for i in range(len(power_set)-1) :
            if power_set[i] :
                power_set[i+1]+=1
                power_set[0]=power_set[i]-1
                if i : power_set[i]=0
                return 1
        return 0

    def spatial_term(powers) :
        """ Returns the properly formatted term corresponding to the
        given powers of xi and eta. """

        if term_type=='arr' :
            if powers[0]==0 : 
                return eta if powers[1]==1 else pow(eta, powers[1])
            if powers[1]==0 : 
                return xi if powers[0]==1 else pow(xi, powers[0])
            return (xi if powers[0]==1 else pow(xi, powers[0]))*(
                eta if powers[1]==1 else pow(eta, powers[1]))
        elif term_type=='coef' :
            if powers[0]==0 : 
                return 'eta' if powers[1]==1 else 'eta%d'%powers[1]
            if powers[1]==0 : 
                return 'xi' if powers[0]==1 else 'xi%d'%powers[0]
            return ('xi' if powers[0]==1 else 'xi'+str(powers[0]))+(
                'eta' if powers[1]==1 else 'eta'+str(powers[1]))
        elif term_type=='expr' :
            if powers[0]==0 : return 'eta^%d'%powers[1]
            if powers[1]==0 : return 'xi^%d'%powers[0]
            return 'xi^%d*eta^%d'%tuple(powers)
        else : raise Error.Input("Unknown term type '%s' in "
                                 "MagnitudeFit.__linear_terms"%term_type)

    def var_term(powers) :
        """ Returns the properly formatted term corresponding to the
        given power(s) of the non-spatial variable(s). """

        if term_type=='arr' :
            part=1.0
            for v,p in zip(var, powers) : 
                if p==1 : part*=v
                elif p!=0 : part*=pow(v,p)
            return part
        else :
            min_i=0
            while powers[min_i]==0 : min_i+=1
            result=var[min_i]
            if powers[min_i]!=1 :
                result+=('%d' if term_type=='coef' else '^%d')%\
                        powers[min_i]
            for i in range(min_i+1,len(var)) :
                if powers[i]!=0 : 
                    if term_type=='expr' : result+='*' 
                    result+=var[i]
                    if powers[i]!=1 :
                        if term_type=='expr' : result+='^'
                        result+=str(powers[i])
            return result

    terms=[]
    for var_order in range(1, max_var_order+1) :
        var_power_set=[(0 if i else var_order) 
                       for i in range(len(var))]
        more_var=True
        while more_var :
            var_part=var_term(var_power_set)
            terms.append(var_part)
            for spatial_order in range(1, max_spatial_order+1) :
                spatial_power_set=[spatial_order, 0]
                more_spatial=True
                while more_spatial :
                    spatial_part=spatial_term(spatial_power_set)
                    if term_type=='arr' :
                        terms.append(var_part*spatial_part)
                    elif term_type=='expr' :
                        terms.append(var_part+'*'+spatial_part)
                    else : terms.append(var_part+spatial_part)
                    more_spatial=next_power_set(spatial_power_set)
            more_var=next_power_set(var_power_set)
    return terms

def all_linear_terms(fit_config, phot=None, subpix_ref=True, 
                     skip_indices=False) :
    """ Returns either a list of numpy arrays containing the combinations
    of independent variables that participate in the fit (if phot is not
    None) or a list of table column names for storing the best fit
    coefficients. """

    if phot is None :
        result=['offset']
    else : 
        result=[ones(len(phot['source'])-(len(skip_indices) if
                                          skip_indices else 0))]
        xi=skip_array(phot['xi'], skip_indices)
        eta=skip_array(phot['eta'], skip_indices)
    if ('spatial' in fit_config.__dict__ and fit_config.spatial) :
        if phot is None : result.extend(linear_terms_subset(
            ['xi', 'eta'], 'xi', 'eta', fit_config.spatial, 0, 'coef'))
        else : result.extend(linear_terms_subset(
            [xi, eta], xi, eta, fit_config.spatial, 0, 'arr'))
    if ('JmK' in fit_config.__dict__ and fit_config.JmK) :
        if phot is None : result.extend(linear_terms_subset(
            ['JmK'], 'xi', 'eta', fit_config.JmK[0], fit_config.JmK[1],
            'coef'))
        else : result.extend(linear_terms_subset(
            [skip_array(phot['J'], skip_indices)-
             skip_array(phot['K'], skip_indices)], xi, eta, 
            fit_config.JmK[0], fit_config.JmK[1], 'arr'))
    for fil in ['r', 'I', 'J', 'K', 'R', 'V', 'B'] :
       # print phot.keys()
      #  print phot['V']
        if (fil in fit_config.__dict__ and
            fit_config.__dict__[fil]) :
            if phot is None : result.extend(linear_terms_subset(
                [fil], 'xi', 'eta', fit_config.__dict__[fil][0],
                fit_config.__dict__[fil][1], 'coef'))
            else : result.extend(linear_terms_subset(
                [skip_array(phot[fil], skip_indices)], xi, eta, 
                fit_config.__dict__[fil][0],
                fit_config.__dict__[fil][1], 'arr'))
    if ('subpix' in fit_config.__dict__ and fit_config.subpix) :
        if phot is None :
            result.extend(linear_terms_subset(
                ['fxI', 'fyI'], 'xi', 'eta', fit_config.subpix[0],
                fit_config.subpix[1], 'coef'))
            last_apply_coef=len(result)
            if subpix_ref :
                result.extend(linear_terms_subset(
                    ['fxR', 'fyR'], 'xi', 'eta', fit_config.subpix[0],
                    fit_config.subpix[1], 'coef'))
        else :
            result.extend(linear_terms_subset(
                [skip_array(phot['x'], skip_indices)%1, 
                 skip_array(phot['y'], skip_indices)%1], xi, eta,
                fit_config.subpix[0], fit_config.subpix[1], 'arr'))
            last_apply_coef=len(result)
            if subpix_ref and ('x ref' in phot) and ('y ref' in phot) :
                result.extend(linear_terms_subset(
                    [skip_array(phot['x ref'], skip_indices)%1, 
                     skip_array(phot['y ref'], skip_indices)%1], xi,
                    eta, fit_config.subpix[0], fit_config.subpix[1], 'arr'))
    else : last_apply_coef=None
    return result, last_apply_coef

class MagnitudeFit :
    """ A base class for all classes doing magnitude fitting against a
    photometric reference (single or master), optionally adding the fitted
    magnitudes as extra columns to the files and optionally writing them to a
    file descriptor. """

    def _update_db(self, values, apind, fit_res, start_src_count,
                   final_src_count) :
        """ Inserts the given values in the database. """

        if self.db is None : return
        args=(self._header['STID'], self._header['FNUM'],
              self._header['CMPOS'], self._header['PROJID'],
              self._header['SPRID'], apind, self._config.version,
              start_src_count, final_src_count, fit_res)+tuple(values)
        statement=('REPLACE INTO `'+self._config.dest_table+
                   '` (`station_id`, `fnum`, `cmpos`, `project_id`, '
                   '`sphotref_id`, `aperture`, `magfit_version`, '
                   '`input_src`, `non_rej_src`, `rms_residuals`, `'+
                   '`, `'.join(self._db_columns)+'`) VALUES (%s'+
                   ', %s'*(len(args)-1)+')')
        sys.stderr.write('statement='+str(statement)+'\n')
        sys.stderr.write('args='+str(args)+'\n')
        self.db(statement, args)

    def _solved(self, num_apertures) :
        """ Returns the best fit coefficients if the solution for the current
        header is already in the database. """

        if self.db is None : return False
        coefficients=[]
        for apind in range(num_apertures) :
            statement=(' FROM `'+self._config.dest_table+
                       '` WHERE `station_id`=%s AND `fnum`=%s AND '
                       '`cmpos`=%s AND `aperture`=%s AND '
                       '`sphotref_id`=%s AND `magfit_version`=%s')
            args=(self._header['STID'], self._header['FNUM'],
                  self._header['CMPOS'], apind, self._header['SPRID'],
                  self._config.version)
            if self.__rederive_fit :
                statement='DELETE '+statement
                self.db(statement, args)
            else :
                statement=('SELECT `'+'`, `'.join(self._db_columns)+'`'+
                           statement)
                ap_coef=self.db(statement, args)
                if ap_coef is None : return False
                else : coefficients.append(None if ap_coef[0] is None else
                                           ap_coef)
        if self.__rederive_fit : return False
        else : return coefficients

    def add_catalogue_info(self, phot) :
        """ Adds the information about the sources in phot that is available
        in the input master catalogue. """
        
        i=0
        new_columns=dict()
        default_cat=(self.__catalogue['default'] 
                     if 'default' in self.__catalogue else None)
        self._default_cat_indices=[]
        while i<len(phot['source']) :
            src=(phot['field'][i], phot['source'][i])
            if src in self.__catalogue : src_dict=self.__catalogue[src]
            elif default_cat : 
                src_dict=default_cat
                self._default_cat_indices.append(i)
            else :
                remove_source_ind(phot, i)
                continue
            for k, v in src_dict.iteritems() :
                if i==0 : new_columns[k]=[]
                new_columns[k].append(v)
            i+=1
        phot.update(new_columns)

    def __format_fiphot_reference(self, ref_phot) :
        """ Formats the given single photometry reference to
        self.__reference: a dictionary indexed by source ids and with values
        a list of magnitudes - one for each aperture. """

        self.__reference=dict()
        for i in range(len(ref_phot['source'])) :
            mag_info=[(ap['mag'][i], ap['mag err'][i]) for ap in 
                      ref_phot['per aperture']]
            bad_source=False
            for mi, ap in zip(mag_info, ref_phot['per aperture']) : 
                bad_source=(bad_source or isnan(mi[0]) or isnan(mi[1]) or
                            mi[1]>self._config.max_mag_err or
                            ap['status flag'][i]!=0)
            if bad_source : continue
            self.__reference[(ref_phot['field'][i], ref_phot['source'][i])]=(
                ref_phot['x'][i], ref_phot['y'][i], mag_info)

    def _match_to_reference(self, phot) :
        """ Returns a photometry structure like phot but for each aperture
        'ref mag' is added - the magnitude the source has in the reference -
        and sources that should not be used in the fit because they are not
        in the reference or do not satisfy the criteria based on catalogue
        quantities are removed. """

        def skip_src(src_ind) :
            """ Returns True iff the source with the given index fails at
            least one catalogue based requirement. """

            if src_ind in no_catalogue : return True
            JmK=phot['J'][src_ind]-phot['K'][src_ind]
            mag=phot[self._config.filchar][src_ind]
            return (mag<self._config.bright_mag or
                    mag>self._config.faint_mag or
                    JmK<self._config.min_JmK or 
                    JmK>self._config.max_JmK or
                    (self._config.AAAonly and 
                     phot['qlt'][src_ind].strip()!='AAA'))
        
        num_ap=len(phot['per aperture'])
        key_list=phot.keys()
        if self._config.reference_subpix : 
            key_list.extend(['x ref', 'y ref'])
        result=dict(zip(key_list, [[] for k in key_list]))
        perapkeys=phot['per aperture'][0].keys()+['ref mag', 'ref mag err']
        for apind in range(num_ap) :
            result['per aperture'].append(dict(zip(perapkeys, [[] for k in
                                                               perapkeys])))
        no_catalogue=set(self._default_cat_indices)
        num_skipped, num_not_in_ref=0, 0
        for src_ind in range(len(phot['source'])) :
            src=(phot['field'][src_ind], phot['source'][src_ind])
            if skip_src(src_ind) : 
                num_skipped+=1
                skipped_example=src
                continue
            if src not in self.__reference :
                num_not_in_ref+=1
                not_in_ref_example=src
                continue
            ref_info=self.__reference[src]
            if self._config.reference_subpix :
                result['x ref'].append(ref_info[0])
                result['y ref'].append(ref_info[1])
                ref_info=ref_info[2]
            for k in phot.keys() :
#                print k
                if k=='per aperture' :
                    for apind, phot_perap in enumerate(phot[k]) : 
                        for apkey, apval in phot_perap.iteritems() :
                            result[k][apind][apkey].append(apval[src_ind])
                        result[k][apind]['ref mag'].append(
                            ref_info[apind][0])
                        result[k][apind]['ref mag err'].append(
                            ref_info[apind][1])
                elif k!='serial' :
                    #print phot[k]
                    result[k].append(phot[k][src_ind])
        if len(result['source']) == 0:
            print "OMGGGGGG"
            print self._header
        return result

    def _get_phot(self) :
        """ Returns the photometry file corresponding to the given header.
        """

        if self.db is None : 
            self._fit_file=self._header['filename']
            f=open('read_files.lst', 'a')
            f.write('reading '+self._header['filename']+'\n')
            f.close()
            return read_fiphot(self._header['filename'])
        missing_keys=dict()
        for k in self.__phot_keys :
            if k not in self._header : missing_keys[k]=self.__phot_keys[k]
        if missing_keys :
            fits=self.db.red_frame(self._header['STID'],
                                   get_night_str(self._header['NIGHT']),
                                   self._header['FNUM'],
                                   self._header['CMPOS'])
            self._header.update(get_keywords(fits, missing_keys))
        self._fit_file=(self.__phot_template.substitute(self._header)+
                        self.db.extensions['rawphot']['fiphot'])
        return read_fiphot(self._fit_file)

    def _append_fitted(self, fitted, src_count) :
        """ Appends the given fitted magnitudes to the input photometry file,
        creates a symlink to the file with the new extension specified at
        construction and updates the calib_status of the file. """

        ap_list=range(len(fitted))
        for f in fitted :
            if f is not None and len(f)!=src_count :
                raise Error.Input(('The number of fitted magnitudes (%d) is '
                                   'not the same as the number of sources '
                                   '(%d) in MagnitudeFit._append_fitted. '
                                   'Fitted magnitudes: %s')%(
                                      len(f), src_count, f))
        name_fmt=self._config.column_name+'[%d]'
        new_columns=[(name_fmt%ap, 
                      ([nan for i in range(src_count)] if fitted[ap] is None 
                       else fitted[ap]), self._config.column_precision) 
                     for ap in ap_list]
        add_columns(self._fit_file, new_columns, self._config.first_column)
        check_photometry(self._fit_file, True, 
                         self._config.check_mfit_columns)
        link_name=(os.path.splitext(self._fit_file)[0]+
                   self._config.file_extension)
        if os.path.lexists(link_name) :
            if os.path.islink(link_name) : os.remove(link_name)
            else : raise Error.File("File '%s' exists and is not a link."%
                                    link_name)
        os.symlink(self._fit_file, link_name)
        if self.db is None : return
        self.db('UPDATE `'+raw_db_table(self._header['IMAGETYP'])+'` SET '
                '`calib_status`=%s WHERE `station_id`=%s AND `fnum`=%s AND '
                '`cmpos`=%s', (self._config.calib_status,
                               self._header['STID'], self._header['FNUM'],
                               self._header['CMPOS']))

    def _downgrade_calib_status(self) :
        """ Decrements the calibration status of the given file to astrometry
        and deletes the raw photometry file. """

        sys.stderr.write('bad photometry encountered:'+str(self._header)+
                         '\n')
        if self.__logconf is not None :
                self.__loglock.acquire()
                self.__logger.critical("Downgrading status of header:"+
                                        str(self._header),
                                        extra=dict(filever=SVNRevision))
                self.__loglock.release()
        sys.stderr.write('downgrading calibration status of'+
                         str(self._header)+'\n')
#        print self._header
        self.db('UPDATE `'+raw_db_table(self._header['IMAGETYP'])+'` SET '
                '`calib_status`=%s WHERE `station_id`=%s AND `fnum`=%s AND '
                '`cmpos`=%s', (object_status['good_astrometry'],
                               self._header['STID'], self._header['FNUM'],
                               self._header['CMPOS']))
        sys.stderr.write('removing:'+self._fit_file+'\n')
        os.remove(self._fit_file)

    def __do_output(self, hat_ids, fitted, formal_errors, phot_flags) :
        """ Writes the appropriate output to the output stream specified at
        construction. """
        return
        
        num_aps=len(fitted)
        assert len(formal_errors)==num_aps
        src_count=len(hat_ids)
        for apind in range(num_aps) : 
            assert(fitted[apind] is None or len(fitted[apind])==src_count)
            assert(len(formal_errors[apind])==src_count)
        fmt=hat_id_fmt+(' %9.5f')*(2*num_aps)
        self.__output_lock.acquire()
        for i in range(src_count) :
            print_args=(hat_ids[i]+
                        tuple((fitted[ap][i] if (fitted[ap] is not None and 
                                                 phot_flags[ap][i]==0) else
                               nan) for ap in range(num_aps))+
                        tuple(formal_errors[ap][i] for ap in range(num_aps)))
            if not reduce(lambda a,b: a or isnan(b), print_args, False) :
                print fmt%print_args
        self.__output_lock.release()

    def __init__(self, reference, config, master_catalogue, db, outlock, 
                 rederive_existing=False, logconf=None, loglock=None) :
        """ Initializes a magnditude fitting thread. 
        
        The input is:
            * reference - the reference against which fitting is done. Should
                          be formatted either as the output of read_fiphot or
                          as a master reference. 
            * master_catalogue - should be a dictionary indexed by sources
                                 (field, source number) containing
                                 dictonaries with relevant 2mass information.
            * db - The calibration database variable. 
            * outlock - A lock to use for ensuring only one thread is
                        outputting at a time """

        self.db=db
        self._config=config
        if self.db is not None :
            self.__phot_template, self.__phot_keys=db.get_name_template(
                db.station, 'rawphot', 'filename')
            self._db_columns=all_linear_terms(
                config, subpix_ref=config.reference_subpix)[0]
        if 'STID' in reference :
            self._header=reference
            self.__format_fiphot_reference(self._get_phot())
        else : self.__reference=reference
        self.__catalogue=master_catalogue
        self.__rederive_fit=rederive_existing
        if logconf is not None and loglock is not None :
            self.__logconf=logconf
            self.__loglock=loglock
        else : self.__logconf=None
        self.__output_lock=outlock

    def __call__(self, header) :
        """ Performs the fit for the frame identified by the given header.
        """

        try :
            if self.__logconf is not None :
                logging.config.fileConfig(os.path.expanduser(self.__logconf))
                self.__logger=logging.getLogger("LogCalib.magfit")
            self._header=header
            phot=self._get_phot()
            if 'source' not in phot or len(phot['source'])==0 : 
                self._downgrade_calib_status()
                return
            src_count=len(phot['source'])
            self.add_catalogue_info(phot)
            coefficients=self._solved(len(phot['per aperture']))
            if not coefficients :
                fit_base=self._match_to_reference(phot)
                if len(fit_base['source'])==0: return
                if self._config.fntype=='linear' :
                    coefficients=self._linear_fit(fit_base)
            fitted=self._apply_fit(phot, coefficients)
            self._append_fitted(fitted, src_count)
            self.__do_output(zip(phot['field'], phot['source']), fitted,
                             [ap['mag err'] for ap in phot['per aperture']],
                             [ap['status flag'] for ap in phot['per aperture']])
        except Exception, ex : 
            print_exc()
            if self.__logconf is not None :
                self.__loglock.acquire()
                self.__logger.critical(str(ex)+"\n"+
                                "".join(format_exception(*sys.exc_info()))+
                                "\nBad header:"+str(header),
                                extra=dict(filever=SVNRevision))
                self.__loglock.release()
            raise

class MagnitudeFitSciPy(MagnitudeFit) :
    """ Implements magnitude fitting using leastsq function of
    scipy.optimize (not thread safe, so it is not done in parallel). """

    def __rejected_indices(self, weighted_fit_diff, weights) :
        """ Returns the indices of those entries in fit_errors which are
        larger than self._config.rej_level*rms(fit_errors) ordered
        large to small, and an estimate of the squared fit residual. """

        fit_diff2=pow(weighted_fit_diff/weights, 2)
        if self._config.error_avg=='weightedmean' :
            res2=mean(pow(weighted_fit_diff,2))/mean(pow(weights,2))
        else :
            avg=eval(self._config.error_avg)
            res2=avg(fit_diff2)
        max_diff2=self._config.rej_level**2*res2
        return filter(lambda i: fit_diff2[i]>max_diff2, 
                      range(len(fit_diff2)-1,-1,-1)), res2

    def _linear_fit(self, phot) :
        """ Performs a linear least squares fit to the given photometry (a
        dictionary with all relevant members) for all apertures, returning
        one list per aperture of the of fit coefficients, which are also
        inserted in the database. 
        
        TODO:
        Perhaps it is worth trying to include the error in the image being
        fitted to the error estimate, in addition to the error in the
        reference at some later time. """

        no_fit_ind=no_fit_indices(phot, self._config.max_mag_err)
        result=[]
        for ap_ind, ap_data in enumerate(phot['per aperture']) :
            skip_ind=no_fit_ind[ap_ind]
            if not skip_ind : skip_ind=False
            predictors=all_linear_terms(self._config, phot,
                                        skip_indices=skip_ind)[0]
            weights=1.0/(skip_array(ap_data['mag err'], skip_ind)+
                         self._config.noise_offset)
            derivatives=map(lambda iv: iv*weights, predictors)
            num_free_coef=len(derivatives)
            mag_array=skip_array(ap_data['mag'], skip_ind)
            mag_difference=(skip_array(ap_data['ref mag'], skip_ind)-
                            mag_array)*weights
            error_func=lambda coef: dot(coef, derivatives)-mag_difference
            deriv_func=lambda coef: derivatives
            initial_guess=zeros(num_free_coef)
            for rej_iter in range(self._config.max_rej_iter) :
                if len(mag_difference)<num_free_coef : 
                    coefficients=None
                    fit_res2=None
                    mag_difference=[]
                    break
                coefficients, covariance, info_dict, msg, status=leastsq(
                    error_func, initial_guess, Dfun=deriv_func,
                    col_deriv=1, full_output=1)
                if status not in [1, 2, 3, 4] :
                    raise Error.Numeric("Linear least squares fitting for "
                                        "aperture %d failed for '%s': %s"%
                                        (ap_ind, self._fit_file, msg))
                bad_ind, fit_res2=self.__rejected_indices(info_dict['fvec'],
                                                          weights)
                if not bad_ind : break
                if rej_iter==self._config.max_rej_iter-1 : break
                derivatives=map(lambda d : delete(d, bad_ind), derivatives)
                mag_difference=delete(mag_difference, bad_ind)
                weights=delete(weights, bad_ind)
            if (self.db is not None) and self._config.dest_table :
                db_coef=([None]*len(predictors) if coefficients is None else
                         coefficients)
                self._update_db(db_coef, ap_ind, sqrt(fit_res2),
                                len(predictors[0]), len(mag_difference))
            result.append(coefficients)
        return result

    def _apply_fit(self, phot, coefficients) :
        """ Returns a structure holding the fitted magnitudes derived by
        using the given fit coefficients. """

        fitted=[]
        predictors, self._last_apply_coef=all_linear_terms(
            self._config, phot, subpix_ref=False)
        for ap_ind in range(len(coefficients)) :
            if coefficients[ap_ind] is None : 
                fitted.append(None)
            else : 
                fitted.append(
                    array(phot['per aperture'][ap_ind]['mag'])+
                    dot((coefficients[ap_ind][:self._last_apply_coef] if
                         self._last_apply_coef else coefficients[ap_ind]),
                        predictors))
        return fitted

def parse_linear_fit_param_str(param_str) :
    """ Parses the given string defining what parameters should participate
    in a linear fit returning a structure with the appropriate members. """

    fit_config=dummy_class()
    for param_val in param_str.split(';') :
        param, value=param_val.split(':')
        if param!='spatial' :
            fit_config.__dict__[param.strip()]=map(eval, value.split(','))
        else : fit_config.__dict__[param.strip()]=eval(value)
    return fit_config

def required_sphotref_keys(db=None) :
    """ Returns a list of keys that should be read from the single
    photometric reference that is sufficient to do any magnitude fitting
    related task. """

    def add_key_strings(required_keys, key_strings) :
        """ Adds the keys in a list of key strings to required_keys. """

        if key_strings is not None :
            map(lambda r: required_keys.update(parse_key_string(r[0])), 
                key_strings)


    required_hdr_keys=dict(
        STID=dict(type='NUMERIC'), FNUM=dict(type='NUMERIC'),
        CMPOS=dict(type='NUMERIC'), NIGHT=dict(type='STRING'),
        PROJID=dict(type='NUMERIC'), EXPTIME=dict(type='NUMERIC'),
        IMAGETYP=dict(type='STRING'), OBJECT=dict(type='STRING'), 
        NAXIS1=dict(type='NUMERIC'), NAXIS2=dict(type='NUMERIC'),
        FILTERS=dict(type='STRING'), NRAC=dict(type='NUMERIC', default=''),
        NDECC=dict(type='NUMERIC', default=''),
        CRA=dict(type='NUMERIC', default=''),
        CDEC=dict(type='NUMERIC', default=''),
        MNTSTATE=dict(type='STRING'))
    if(db is not None) :
        for temp_name in ['masterphotref', 'master2MASS', 'sphotrefstat',
                          'mphotrefstat'] :
            required_hdr_keys.update(db.get_name_template(
                db.station, 'magfit', temp_name)[1])
        add_key_strings(required_hdr_keys,
                        db('SELECT `config_keys` FROM `RawPhotCfg`',
                           no_simplify=True))
        add_key_strings(required_hdr_keys,
                        db('SELECT `keys` FROM `SinglePhotRef`',
                           no_simplify=True))
    for k in ['APERTURE', 'PHOT_TYPE'] :
        if k in required_hdr_keys : del required_hdr_keys[k]
    return required_hdr_keys

def read_config_file(fname) :
    """ Reads the configuration for the magnitude fitting from the given 
    file.

    The file should have the following format:
        [section name]
        key1 = value1
        key2 = value2
        ...
        [section name]
        ...
    Where [] is a valid section resulting in the defined keys being directly
    members of the configuration and can be omitted if this is the first
    section in the file. """

    f=open(fname, 'r')
    result=dummy_class()
    destination=result.__dict__
    for l in f :
        clean_line=l.strip()
        if clean_line=='' : continue
        if clean_line[0]=='[' :
            assert(clean_line[-1]==']')
            section=clean_line[1:-1].strip().replace(' ', '_')
            if section in ['magfit.single', 'magfit.master'] :
                if 'magfit' not in result.__dict__ :
                    result.magfit=dummy_class()
                if section[7:] not in result.magfit.__dict__ :
                    result.magfit.__dict__[section[7:]]=dummy_class()
                destination=result.magfit.__dict__[section[7:]].__dict__
            elif section=='' : destination=result.__dict__
            else :
                if section not in result.__dict__ :
                    result.__dict__[section]=dummy_class()
                destination=result.__dict__[section].__dict__
        else :
            key, value=l.split('=',1)
            key=key.strip()
            value=value.strip()
            destination[key]=eval(value)
    result.mphotref.rms_fit_param=parse_linear_fit_param_str(
        result.mphotref.rms_fit_param)
    result.magfit.single.__dict__.update(
        parse_linear_fit_param_str(result.magfit.single.param_str).__dict__)
    result.magfit.master.__dict__.update(
        parse_linear_fit_param_str(result.magfit.master.param_str).__dict__)
    return result

def get_config(sphotref_header, db) :
    """ Returns a configuration structure defining various parameters
    necessary for generating a master photometric reference. """

    def read_mphotref(config) :
        """ Adds to config mebers containing all the configuaration contained
        in the MasterPhotRefCfg table. """
        
        querry=('SELECT `fov_safety_fac`, `fov_alarm`, `rej_outliers`, '
                '`rej_iterations`, `master_cat_cols`, `cat_point_round`, '
                '`min_measurements`, `min_meas_med`, `grcollect_tempdir`, '
                '`rms_fit_param`, `rms_fit_err_avg`, `rms_fit_rej_lvl`, '
                '`max_rms_above_fit`, `max_rms_quantile`, `rms_fit_bright`, '
                '`version` FROM `MasterPhotRefCfg` WHERE `project_id`=%s AND'
                ' `sphotref_id`=%s ORDER BY `version` DESC LIMIT 1')
        record=db(querry, (sphotref_header['PROJID'],
                           sphotref_header['SPRID']))
        if record is None : 
            record=db(querry, (sphotref_header['PROJID'], 0))
        if record is None :
            raise Error.Database('No specific or default MasterPhotRefCfg '
                                 'entry found for project ID %d, sphotref ID'
                                 ' %d.'%(sphotref_header['PROJID'],
                                         sphotref_header['SPRID']))
        (config.mcat.fov_safety_fac, config.mcat.fov_alarm,
         config.mphotref.rej_outliers, config.mphotref.rej_iterations, 
         master_cat_columns, config.mcat.round_pointing,
         config.mphotref.min_measurements, config.mphotref.min_meas_med,
         config.mphotref.grcollect_tempdir, rms_fit_param_str,
         rms_fit_err_avg, config.mphotref.rms_fit_rej_lvl,
         config.mphotref.max_rms_above_fit, config.mphotref.max_rms_quantile,
         config.mphotref.rms_fit_bright_mag_min,
         config.mphotref.version)=record
        config.mphotref.rms_fit_err_avg=eval(rms_fit_err_avg)
        config.mcat.columns=master_cat_columns.strip().split(',')
        config.mphotref.rms_fit_param=parse_linear_fit_param_str(
            rms_fit_param_str)

    def read_rawphot(config) :
        """ Appends the configuration from the RawPhotCfg table to the given
        configuration structure. """
        
        records=db('SELECT `config_cond`, `bright_mag`, `faint_mag`, '
                   '`apertures`, `version` FROM `RawPhotCfg` WHERE '
                   '`project_id`=%s AND `station_id`=%s and `version`=%s',
                   (sphotref_header['PROJID'], sphotref_header['STID'],
                    raw_phot_cfg_version(sphotref_header['STID'],
                                         sphotref_header['PROJID'], db)), 
                   no_simplify=True)
        for r in records :
            if eval(Template(r[0]).substitute(sphotref_header)) :
                (config.mcat.bright_mag, config.mcat.faint_mag, aperture_str,
                 config.rawphot_ver)=r[1:]
                config.num_apertures=aperture_str.count(',')+1
                return

    def read_magfit(config) :
        """ Returns a class containing the configurations of how to perform
        magnitude fitting for a single and master modes. """

        cfg=dummy_class()
        for mode in ['single', 'master'] :
            cfg.__dict__[mode]=dummy_class()
            subcfg=cfg.__dict__[mode]
            querry=('SELECT `func_type`, `func_param`, `AAAonly`, '
                    '`bright_mag_min`, `min_JmK`, `max_JmK`, '
                    '`ref_frame_stars`, `noise_offset`, `max_mag_err`, '
                    '`rej_level`, `max_rej_iter`, `error_avg`, '
                    '`count_weight`, `version` FROM `MagFitCfg` WHERE '
                    '`mode`=%s AND `project_id`=%s AND `sphotref_id`=%s '
                    'ORDER BY `version` DESC LIMIT 1')
            record=db(querry, (mode, sphotref_header['PROJID'],
                               sphotref_header['SPRID']))
            if record is None :
                record=db(querry, (mode, sphotref_header['PROJID'], 0))
            if record is None :
                raise Error.Database('No specific or default MagFitCfg '
                                     'entry found for project ID %d, '
                                     'sphotref ID %d.'%(
                                         sphotref_header['PROJID'],
                                         sphotref_header['SPRID']))
            (subcfg.fntype, param_str, subcfg.AAAonly, subcfg.bright_mag_min,
             subcfg.min_JmK, subcfg.max_JmK, subcfg.ref_frame_stars,
             subcfg.noise_offset, subcfg.max_mag_err, subcfg.rej_level,
             subcfg.max_rej_iter, subcfg.error_avg, subcfg.count_weight,
             subcfg.version)=record
            subcfg.__dict__.update(
                parse_linear_fit_param_str(param_str).__dict__)
            subcfg.column_name=mode[0]+'prmag'
            subcfg.column_precision=1e-5
            subcfg.file_extension=db.extensions['magfit'][mode+'photref']
            subcfg.calib_status=object_status[mode+'_magfit']
            subcfg.dest_table=mode.capitalize()+'PhotRefMagFit'
            subcfg.reference_subpix=(mode=='single')
        cfg.single.first_column=0
        cfg.master.first_column=config.num_apertures
        config.magfit=cfg

    def add_name_templates(config) :
        """ Adds the appropriate options to the given configuration structure
        defining the various names for files and the keys that must be read
        from the single photometry header needed to substitute in those
        names. """

        def read_template(temp_name) :
            """ Reads the temp_name template from the magfit component and
            returns it. """

            return db.get_name_template(db.station, 'magfit', temp_name)[0]
        
        config.mphotref.template=read_template('masterphotref')
        config.mcat.template=read_template('master2MASS')
        config.magfit.single.stat_template=read_template('sphotrefstat')
        config.magfit.master.stat_template=read_template('mphotrefstat')

    config=dummy_class()
    config.mphotref=dummy_class()
    config.mcat=dummy_class()
    config.magfit=dummy_class()
    read_mphotref(config)
    read_rawphot(config)
    read_magfit(config)
    add_name_templates(config)
    return config

def read_master_catalogue(fname, columns, str_id=False) :
    """ Returns a source indexed dictionary of all sources in the given 2MASS
    querry. """

    input_2mass=open(fname, 'r')
    id_col=columns.index('id')
    result=dict()
    for l in input_2mass :
        if l.startswith('#') : continue
        entries=l.split()
        src_id=(entries[id_col] if str_id else parse_hat_id(entries[id_col]))
        result_entry=dict()
        for i,k in enumerate(columns) :
            if i==id_col : continue
            elif k=='qlt' : v=entries[i]
            else : v=eval(entries[i])
            result_entry[k]=v
        result[src_id]=result_entry
    input_2mass.close()
    return result

class SinglePhotRef :
    """ A callable class that returns the relevant header keywords from the
    single photometric reference that should be used for magnitude fitting of
    a frame with a given header. """

    def __init__(self, db, project, field, cmpos, imagetype, extra_hdr_keys):
        """ Prepare as much as possible to return quickly the suitable single
        photometric reference for frames with the given restrictions. """

        frame_table=raw_db_table(imagetype)
        candidates=db('SELECT `S`.`sphotref_id`, `F`.`station_id`, '
                      '`F`.`night`, `F`.`fnum`, `S`.`condition` FROM '
                      '`SinglePhotRef` AS `S` JOIN `'+frame_table+
                      '` AS `F` ON (`S`.`station_id`=`F`.`station_id`'
                      ' AND `S`.`fnum`=`F`.`fnum` AND '
                      '`S`.`cmpos`=`F`.`cmpos`) WHERE '
                      '`S`.`project_id`=%s AND `S`.`object`=%s AND '
                      '`S`.`cmpos`=%s AND `S`.`enabled`=%s', 
                      (project, field, cmpos, 1), no_simplify=True)
        need_keys=db.get_name_template(db.station, 'astrometry','outpath')[1]
        phot_keys=db.get_name_template(db.station, 'rawphot', 'filename')[1]
        need_keys.update(phot_keys)
        need_keys.update(extra_hdr_keys)
        for k in ['STID', 'NIGHT', 'FNUM', 'CMPOS', 'IMAGETYP', 'OBJECT',
                  'PROJID', 'SPRID'] : 
            if k in need_keys : del need_keys[k]
        self.__headers, self.__conditions=[], []
        for photref_id, stid, night, fnum, condition in candidates :
            night_str=get_night_str(night)
            header=dict(STID=stid, NIGHT=night_str, FNUM=fnum, CMPOS=cmpos,
                        IMAGETYP=imagetype, OBJECT=field, PROJID=project,
                        SPRID=photref_id)
            fits=db.red_frame(stid, night_str, fnum, cmpos)
            if need_keys : header.update(get_keywords(fits, need_keys))
            patch_header_keywords(header)
            self.__headers.append(header)
            self.__conditions.append(Template(condition))

    def __call__(self, frame_header_or_id) :
        """ Returns the photometric reference that is suitable for the given
        header (if frame_header_or_id is a dictionary of header keywords) or
        it has the given id if it is an integer. """

        assert (type(frame_header_or_id) is dict or
                type(frame_header_or_id) is int) 
        found=None
        for condition, header in zip(self.__conditions, self.__headers) :
            if type(frame_header_or_id) is dict : 
                frame_header_or_id['SPRID']=header['SPRID']
            if ((type(frame_header_or_id) is int and
                 header['SPRID']==frame_header_or_id) or 
                (type(frame_header_or_id) is dict and
                 eval(condition.substitute(frame_header_or_id)))) :
                if found is not None : 
                    raise Error.Database('Multiple single photometric'
                                         'references found for header '+
                                         str(frame_header_or_id)+
                                         '%d and %d'%(found['SPRID'], 
                                                      header['SPRID']))

                found=header
        return found

def get_work(db, project, sphotref_id, field, cmpos, imagetype, mode) :
    """ Returns a list of the frames for which magnitude fitting is to be
    performed. """

    def sphotref_condition() :
        """ Returns the condition that the headers of frames should satisfy
        in order to be processed in this run. """

        return Template(db('SELECT `condition` FROM `SinglePhotRef` WHERE '
                           '`sphotref_id`=%s AND `project_id`=%s',
                           (sphotref_id, project))[0])

    def required_keys() :
        """ A dictionary of keys that must be read from raw frame headers in
        order to test the single photometric reference condition. """

        keys=parse_key_string(db('SELECT `keys` FROM `SinglePhotRef` WHERE'
                                 ' `sphotref_id`=%s AND `project_id`=%s',
                                 (sphotref_id, project))[0])
        for k in ['PROJID', 'OBJECT', 'CMPOS', 'STID', 'FNUM', 'NIGHT',
                  'IMAGETYP', 'SPRID'] :
            if k in keys : del keys[k]
        return keys

    def get_hdr(record, extra_keys) :
        """ Returns a header like structure generated from the entries of the
        given record. """

        hdr=dict(PROJID=record[0], OBJECT=record[1], CMPOS=record[2],
                 STID=record[3], FNUM=record[4],
                 NIGHT=get_night_str(record[5]), IMAGETYP=imagetype,
                 SPRID=sphotref_id)
        fname=db.red_frame(stid=hdr['STID'], fnum=hdr['FNUM'],
                           cmpos=hdr['CMPOS'], night=hdr['NIGHT'])
        hdr.update(get_keywords(fname, extra_keys))
        return hdr

    frame_tbl=raw_db_table(imagetype)
    condition=sphotref_condition()
    keys=required_keys()
    querry_start=('SELECT `project_id`, `object`, `cmpos`, `station_id`, '
                  '`fnum`, `night` FROM `'+frame_tbl+'` WHERE '
                  '`project_id`=%s AND `object`=%s AND `cmpos`=%s')
    querry_end='ORDER BY `station_id`, `fnum`'
    frames_todo=db(querry_start+' AND `calib_status`>=%s AND '
                   '`calib_status`<%s '+querry_end, 
                   (project, field, cmpos, object_status['raw_photometry'],
                    object_status[mode+'_magfit']), no_simplify=True)
    if frames_todo is None : frames_todo=[]
    else : frames_todo=map(lambda f: get_hdr(f, keys), frames_todo)
    frames_todo=filter(lambda f: eval(condition.substitute(f)), frames_todo)
    frames_done=db(querry_start+' AND `calib_status`=%s '+querry_end, 
                   (project, field, cmpos, object_status[mode+'_magfit']),
                   no_simplify=True)
    if frames_done is None : frames_done=[]
    else : frames_done=map(lambda f: get_hdr(f, keys), frames_done)
    frames_done=filter(lambda f: eval(condition.substitute(f)), frames_done)
    return frames_todo, frames_done

def read_work_list(filename, required_header_keys) :
    """ Returns a list of headers with a 'filename' keyword for the frames that
    listed in the given file. """

    f=open(filename, 'r')
    result=[]
    for l in f :
        phot, fits=l.split()
        header=get_image_vars(fits, required_header_keys)
        header['filename']=phot
        result.append(header)
    f.close()
    return result

def start_grcollect(mode, config, sph_header, input_stream) :
    """ Returns the grcollect command to use for generating statistics from
    the magnitude fit. """

    stat_template=(config.magfit.single.stat_template if mode=='single' else 
                   config.magfit.master.stat_template)
    magfit_stat=stat_template.substitute(sph_header)
    grcollect_cmd=['grcollect', '-', '-V', '--stat']
    stat_columns=range(2,2*config.num_apertures+2)
    if config.mphotref.rej_outliers :
        out_columns=['count', 'rcount', 'rmedian', 'rmediandev', 
                     'rmedianmeddev']
        grcollect_cmd.append(','.join(out_columns))
        for col in stat_columns :
            grcollect_cmd.extend(
                ['--rejection', 'column=%d,iterations=%d,median,%s'%
                 (col, config.mphotref.rej_iterations,
                  config.mphotref.rej_outliers)])
    else : 
        out_columns=['count', 'median', 'mediandev', 'medianmeddev']
        grcollect_cmd.append('count,count,median,mediandev,medianmeddev')
    out_columns=['id']+out_columns
    grcollect_cmd.extend(['--col-base', '1', '--col-stat', 
                          ','.join(map(str, stat_columns)),
                          '--max-memory', '2g', '--tmpdir',
                          config.mphotref.grcollect_tempdir, '--output',
                          magfit_stat])
    if not os.path.exists(config.mphotref.grcollect_tempdir) :
        os.makedirs(config.mphotref.grcollect_tempdir)
    return Popen(grcollect_cmd, stdin=input_stream, stderr=PIPE), out_columns

def init_magfit(mode, sph_header, db, config, logconf) :

    def parse_photref_file(f) :
        """ Reads the relevant columns from the given master photometry
        reference file and returns a dictionary indexed by source of
        (magnitude, magnitude scatter). """

        ref=dict()
        for line in f :
            if line[0]=='#' : continue
            entries=line.split()
            src_id=parse_hat_id(entries[0])
            ref[src_id]=tuple(map(eval, entries[3:5]))
        return ref

    def get_master_photref() :
        """ Returns the master photometry reference that should be used for
        the given magnitude fitting task. """

        reference=dict()
        template_keys=dict(sph_header)
        ap_indices=range(config.num_apertures)
        ap_ref=[None for ap in ap_indices]
        all_sources=[]
        for ap in ap_indices :
            template_keys['APERTURE']=ap
            mphot=open(config.mphotref.template.substitute(
                template_keys), 'r')
            ap_ref[ap]=parse_photref_file(mphot)
            mphot.close()
            all_sources.extend(ap_ref[ap].keys())
        all_sources=sorted(list(set(all_sources)))
        for src in all_sources :
            reference[src]=[ap_ref[ap][src] if src in ap_ref[ap] else 
                            (nan, nan) for ap in ap_indices]
        return reference

    def faint_mag(mfit_cfg, master_cat, reference=None) :
        """ Returns the faint magnitude limit which will result in at most
        config.magfit.ref_frame_stars count of sources that fall within the
        single photometry reference with the selection. 
        
        If reference is given, only sources from the reference are counted
        toward the limit. """

        def good(src) :
            """ Returns True iff the given source passes all catalogue based
            requirements and is not missing from the reference. """

            if reference is not None and src not in reference : return False
            cat_info=master_cat[src]
            JmK=cat_info['J']-cat_info['K']
            filchar=mfit_cfg.filchar
            return (cat_info[filchar]>mfit_cfg.bright_mag and
                    JmK>mfit_cfg.min_JmK and JmK<mfit_cfg.max_JmK and
                    (not mfit_cfg.AAAonly or cat_info['qlt'].strip()=='AAA'))

        fchar=mfit_cfg.filchar
        mags=map(lambda s: master_cat[s][fchar], 
                 filter(good, master_cat.keys()))
        if len(mags)<=mfit_cfg.ref_frame_stars : return inf
        else : 
            mags.sort()
            return (mags[mfit_cfg.ref_frame_stars-1]+
                    mags[mfit_cfg.ref_frame_stars])/2.0

    master_cat=read_master_catalogue(
        config.mcat.template.substitute(sph_header), config.mcat.columns)
    fit_config=config.magfit.__dict__[mode]
    fit_config.filchar=sph_header['FILTERS']
    fit_config.bright_mag=mag_from_mag_per_minute(fit_config.bright_mag_min, 
                                                  sph_header['EXPTIME'])
    if db is None :
        phot_fname=sph_header['filename']
    else : phot_fname=(db.get_name_template(
        db.station, 'rawphot', 'filename')[0].substitute(sph_header)+
                db.extensions['rawphot']['fiphot'])
    sph_phot=read_fiphot(phot_fname, common_only=True)
    sph_phot_stars=set(zip(sph_phot['field'], sph_phot['source']))
    if mode=='single' :
        reference=sph_header
        fit_config.faint_mag=faint_mag(fit_config, master_cat,sph_phot_stars)
        fit_config.check_mfit_columns=False
    else :
        reference=get_master_photref()
        ref_stars=set(reference.keys()) & sph_phot_stars
        fit_config.faint_mag=faint_mag(fit_config, master_cat, ref_stars)
        fit_config.check_mfit_columns=True
    master_cat['default']=dict.fromkeys(config.mcat.columns, nan)
    del master_cat['default']['id']
    manager=Manager()
    return MagnitudeFitSciPy(reference=reference, config=fit_config,
                             master_catalogue=master_cat, db=db,
                             outlock=manager.Lock(), logconf=logconf, 
                             loglock=manager.Lock())

def mag_from_mag_per_minute(mag_min, exptime) :
    """ Returns the magnitude scaled by the given exposure time corresponding
    to the given magnitude in a 1 min exposure. """

    return mag_min+2.5*log10(exptime/60.0)

def start_magfit(project, sphotref_id, field, cmpos, imagetype, options, db,
                 mode) :
    """ Starts magnitude fitting for the given field/camera/imagetype using
    the number of processes specified in options. """

    command=['MagnitudeFitting.py', mode, str(project), str(sphotref_id),
             field, str(cmpos), imagetype, '--database', 
             os.path.join(os.path.expanduser(options.config_dir),
                          options.localdb),
             '--processes', str(options.max_threads)]
    magfit_ver=db('SELECT MAX(`version`) FROM `MagFitCfg` WHERE '
                  '`project_id`=%s AND `sphotref_id`=%s AND `mode`=%s',
                  (project, sphotref_id, mode))[0]
    if magfit_ver is None :
        magfit_ver=-db('SELECT MAX(`version`) FROM `MagFitCfg` WHERE '
                       '`project_id`=%s AND `sphotref_id`=%s AND `mode`=%s',
                       (project, 0, mode))[0]
    return Popen(command, stdout=PIPE, stderr=PIPE), magfit_ver

def log_magfit(db, imagetype, field, cmpos, logger, 
               mode) :
    """ Generates appropriate log messages about how the magnitude fitting
    went, based on the database. """

    frame_tbl=raw_db_table(imagetype)
    magfit_tbl=mode.capitalize()+'PhotRefMagFit'
    statistics=db('SELECT `F`.`aperture`, `F`.`station_id` IS NULL, '
                  '`F`.`non_rej_src`="0", COUNT(*) FROM `'+frame_tbl+
                  '` AS `I` LEFT JOIN `'+magfit_tbl+'` AS `F`'
                  ' ON (`I`.`station_id`=`F`.`station_id` AND '
                  '`I`.`fnum`=`F`.`fnum` AND `I`.`cmpos`=`F`.`cmpos`) '
                  'WHERE `I`.`object`=%s AND `I`.`cmpos`=%s AND '
                  '`I`.`calib_status`>=%s GROUP BY `F`.`aperture`, '
                  '`F`.`station_id` IS NULL, `F`.`non_rej_src`="0"', 
                  (field, cmpos, object_status['raw_photometry']),
                  no_simplify=True)
    current_ap=None
    formatted_stat=dict()
    total_attempted=0
    good_count, bad_count=0,0
    for (ap, not_done, failed, count) in statistics :
        if ap!=current_ap :
            if current_ap is not None :
                total_attempted=max(good_count+bad_count,total_attempted)
            formatted_stat[current_ap]=(good_count, bad_count)
            good_count, bad_count=0,0
            current_ap=ap
        if failed==1 or not_done==1 : bad_count=count
        elif failed==0 : good_count=count
    formatted_stat[current_ap]=(good_count, bad_count)
    good_count, bad_count=0,0
    if None in formatted_stat :
        total_skipped=sum(formatted_stat[None])
        if total_skipped!=0 :
            logger.error('Magnitude fit was not attempted for any aperture '
                         'of %d out of %d %s frames from field %s camera %d'
                         %(total_skipped, total_attempted+total_skipped,
                           imagetype, field, cmpos))
    for ap, (good, bad) in sorted(formatted_stat.items()) :
        if ap is None : continue
        if bad : logger.warning(
            'Magnitude fitting for aperture %d failed for %d out of %d '
            'attempted %s frames from field %s camera %d'%
            (ap, bad, good+bad, imagetype, field, cmpos))

def output_done(hdr_list, mode, db, num_aps, logger=None, destination=None) :
    """ Writes the appropriate magnitudes and errors to the output stream. 
    """

    if db is not None :
        phot_template, phot_keys=db.get_name_template(db.station, 'rawphot',
                                                      'filename')
        missing_keys=dict()
        for k in phot_keys :
            if k not in hdr_list[0] : missing_keys[k]=phot_keys[k]
        file_ext=db.extensions['magfit'][mode+'photref']
    fmt=hat_id_fmt+(' %9.5f')*(2*num_aps)
    fit_col=mode[0]+'prmag[%d]'
    for header in hdr_list :
        if db is None :
            try :
                 phot=read_fiphot(header)
            except :
                 print 'faile to read', header
                 raise
        else :
            if missing_keys :
                fits=db.red_frame(header['STID'],
                                  get_night_str(header['NIGHT']),
                                  header['FNUM'], header['CMPOS'])
                header.update(get_keywords(fits, missing_keys))
            phot=read_fiphot(phot_template.substitute(header)+file_ext)
        src_count=len(phot['source'])
        if(logger is not None) :
            logger.debug("Outputting %d sources.\n"%src_count)
        for i in range(src_count) :
            print_args=((phot['field'][i], phot['source'][i])+
                        tuple(phot[fit_col%ap][i] for ap in range(num_aps))+
                        tuple(phot['per aperture'][ap]['mag err'][i] 
                              for ap in range(num_aps)))
            if not reduce(lambda a,b: a or isnan(b), print_args, False) :
                if destination is None : print fmt%print_args
                else : destination.write(fmt%print_args+'\n')

def patch_header_keywords(header) :
    """ Converts RA from hours to degrees and handles the possibility that
    CRA is in the header instead of NRAC ond CDEC instead of NDECC. """

    if 'NRAC' in header :
        for k in ['NRAC', 'CRA'] : 
            if header[k]!='' : header[k]*=15.0
        if header['NRAC']=='' : 
            if header['CRA']=='' : 
                raise Error.FitsHeader('Neither NRAC nor CRA header '
                                       'keywords found in %s'%fits)
            else : header['NRAC']=header['CRA']
    if 'NDECC' in header :
        if header['NDECC']=='' : 
            if header['CDEC']=='' :
                raise Error.FitsHeader('Neither NDECC nor CDEC '
                                       'header keywords found in %s'%
                                       fits)
            else : header['NDECC']=header['CDEC']

def get_phot_type(phot_list) :
    """ Returns the photometry type for the given file list, verifying
    uniqueness. """
    
    if type(phot_list[0]) is dict :
        phot_filenames=map(lambda d: d['filename'], phot_list)
    else : phot_filenames=phot_list
    phot_types=set(map(
        lambda fname: os.path.basename(os.path.dirname(fname)),
        phot_filenames))
   # assert len(phot_types)==1
    return phot_types.pop()

def parse_command_line() :
    """ Parses the command line options. """

    parser=OptionParser(usage='%prog <single|master> <field> <cmpos> '
                       '<imagetype>',
                        description='Performs magnitude fitting for all '
                        'frames with raw photometry from the given field and'
                        ' camera position.', 
                        version="%prog SVN Revsion: "+str(SVNRevision))
    parser.add_option("--database", action="store", type="string", 
                      default=os.path.expanduser('~/CALIBLOG/.calib_db'),
                      help='A file to read the options for connecting to the'
                      ' database. Default: %default.')
    parser.add_option('-p', '--processes', action='store', type='int',
                      default='6', help='The number of processes to use. '
                      'Default: %default')
    parser.add_option('--log-config', action='store', type='string',
                      dest='log_conf', default='~/CALIBLOG/logging.conf', 
                      help="A file containing the configuration of the "
                      "logger (%default by default). See python logging "
                      "module documentation for the format. The name of the "
                      "logger should be 'LogCalib'.")
    parser.add_option('--manual-frame-list', action='store', default='',
                      help='The list of files to perform magnitude fitting '
                      'from is read from the given file. With this option '
                      'the positional arguments should be <HATNet|HATSouth> '
                      '<single|master> <photometric reference fits file> '
                      '<photometric reference photometry file>. If passed, '
                      'no database updates are done and statistics is '
                      'collected only if the --stat-file option is '
                      'supplied.')
    parser.add_option('--config-file', action='store', default='magfit.cfg',
                      help='The name of a file to read various configuration'
                      ' options from. Only used if --manual-frame-list '
                      'options is passed. Default: %default.')
    options, positional_arguments=parser.parse_args()
    if(options.manual_frame_list) :
        network, options.mode, options.sphotref_fits, options.sphotref_phot=\
                positional_arguments
        if network.upper()=='HATNET' : HATUtil.Stations=HATNetStations
        elif network.upper()=='HATSOUTH' : HATUtil.Stations=HATSouthStations
        else : raise Error.CommandLine("Unrecognized telescope network: "+
                                       network)
    else : 
        (options.mode, options.project, options.sphotref_id, options.field,
         options.cmpos, options.imagetype)=positional_arguments
        options.sphotref_id=eval(options.sphotref_id)
        options.cmpos=eval(options.cmpos)
        options.project=eval(options.project)
    return options

if __name__=='__main__' :
    options=parse_command_line()
    if(not options.manual_frame_list) :
        db=get_db(options.database)
        required_hdr_keys=required_sphotref_keys(db)
        sph_header=SinglePhotRef(db, options.project, options.field,
                                 options.cmpos, options.imagetype,
                                 required_hdr_keys)(options.sphotref_id)
        if sph_header is None : raise Error.Database(
            'No single photometric reference found for project %d with id %d'
            %(project, sphotref_id))
        config=get_config(sph_header, db)
        pending, done=get_work(db, options.project, options.sphotref_id,
                               options.field, options.cmpos,
                               options.imagetype, options.mode)
    else :
        required_hdr_keys=required_sphotref_keys()
        sph_header=get_image_vars(options.sphotref_fits, required_hdr_keys)
        sph_header['filename']=options.sphotref_phot
        patch_header_keywords(sph_header)
        config=read_config_file(options.config_file)
        pending=read_work_list(options.manual_frame_list, required_hdr_keys)
#        sph_header['PHOT_TYPE']=get_phot_type(pending)
        done=[]
        db=None
    if pending :
        magfit=init_magfit(options.mode, sph_header, db, config,
                           options.log_conf)
        workers=Pool(options.processes)
        workers.map_async(
            magfit, pending,
            chunksize=max(10, len(pending)/(10*options.processes))).get()
        workers.close()
        workers.join()

    sys.stderr.write('Outputting %d pre-fitted frames.\n'%len(done))
    logging.config.fileConfig(os.path.expanduser(options.log_conf))
    logger=logging.getLogger("LogCalib.magfit")
    logger.info('Outputting %d pre-fitted frames.\n'%len(done), 
                extra=dict(filever=SVNRevision))
    if done : output_done(done, options.mode, db, config.num_apertures,
                          logger)
