#!/usr/bin/python
# encoding: utf-8

##Calculation function modified based on STRING Database scripts.

###########################################################################################################
## translate to Python from NEWSTRING_compute_scores.pm
##
## This package should serve as the only point where STRING/STITCH 'combined_score' values are computed.
## It should do all the appropriate homology corrections etc the code is specifically put into this
## module to avoid inconsistencies between frontend and backend.
############################################################################################################

## all the routines below are assuming that scores are provided as values between 0 and 1,
## and are returning a combined score between 0 and 1.

import sys
import os


############################################################################################################
## subroutine: compute_combined_score_protein_protein ()
##
############################################################################################################

## same prior for all species. this value is locally stored here,
## but can be overriden if given as the last parameter of the various functions below.

## for STRING version 8.1 / STITCH version 2.0:
local_prior_pp = 0.063 ## note: when Michael recalculated this prior, it went to 0.052. This will, however, not be changed for STRING 8.1
local_prior_pc = 0.029
local_prior_cc = 0.026

def compute_combined_score_protein_protein (
        nscore, nscore_transferred,
        fscore,
        pscore,
        hscore,
        ascore, ascore_transferred,
        escore, escore_transferred,
        dscore, dscore_transferred,
        tscore, tscore_transferred,
        prior = None):

    """
    >>> args = ", ".join( ["0.6"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["1.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["0.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(0.999*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(1.001*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    """

    if prior is None:
        prior = local_prior_pp

    ## some basic checks to learn about impossible errors the hard way:

    assert prior > 0.01, "error: prior prior far too low !!"
    assert prior < 0.50, "error: prior prior far too high !!"

    ## next, remove the priors from all scores. notabene: no prior correction on the hscore is needed.

    nscore_prior_corrected             = compute_prior_away (nscore, prior)
    nscore_transferred_prior_corrected = compute_prior_away (nscore_transferred, prior)
    fscore_prior_corrected             = compute_prior_away (fscore, prior)
    pscore_prior_corrected             = compute_prior_away (pscore, prior)
    ascore_prior_corrected             = compute_prior_away (ascore, prior)
    ascore_transferred_prior_corrected = compute_prior_away (ascore_transferred, prior)
    escore_prior_corrected             = compute_prior_away (escore, prior)
    escore_transferred_prior_corrected = compute_prior_away (escore_transferred, prior)
    dscore_prior_corrected             = compute_prior_away (dscore, prior)
    dscore_transferred_prior_corrected = compute_prior_away (dscore_transferred, prior)
    tscore_prior_corrected             = compute_prior_away (tscore, prior)
    tscore_transferred_prior_corrected = compute_prior_away (tscore_transferred, prior)


    ## then, combine the direct and transferred scores for each category:

    nscore_both_prior_corrected = 1.0 - (1.0 - nscore_prior_corrected) * (1.0 - nscore_transferred_prior_corrected)
    ascore_both_prior_corrected = 1.0 - (1.0 - ascore_prior_corrected) * (1.0 - ascore_transferred_prior_corrected)
    escore_both_prior_corrected = 1.0 - (1.0 - escore_prior_corrected) * (1.0 - escore_transferred_prior_corrected)
    dscore_both_prior_corrected = 1.0 - (1.0 - dscore_prior_corrected) * (1.0 - dscore_transferred_prior_corrected)
    tscore_both_prior_corrected = 1.0 - (1.0 - tscore_prior_corrected) * (1.0 - tscore_transferred_prior_corrected)

    ## now, do the homology correction on pscore and tscore:

    phscore_prior_corrected = pscore_prior_corrected * (1.0 - hscore)
    thscore_prior_corrected = tscore_both_prior_corrected * (1.0 - hscore)

    ## next, do the 1 - multiplication:

    combined_score_one_minus = (
        (1.0 - nscore_both_prior_corrected) *
        (1.0 - fscore_prior_corrected) *
        (1.0 - phscore_prior_corrected) *
        (1.0 - ascore_both_prior_corrected) *
        (1.0 - escore_both_prior_corrected) *
        (1.0 - dscore_both_prior_corrected) *
        (1.0 - thscore_prior_corrected)
    )

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                                 ## scale down
    combined_score += prior                                         ## and add prior.

    ## out-of-bounds check once again, to make sure:
    if combined_score < 0.0: combined_score = 0.0
    if combined_score > 1.0: combined_score = 1.0

    ## and return the result.

    return combined_score

def compute_combined_score_protein_protein_new_full (
        cscore,
        nscore,
        fscore,
        pscore,
        hscore,
        ascore,
        escore,
        dscore,
        tscore,
        prior = None):

    """
    >>> args = ", ".join( ["0.6"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["1.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["0.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(0.999*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(1.001*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    """

    if prior is None:
        prior = local_prior_pp

    ## some basic checks to learn about impossible errors the hard way:

    assert prior > 0.01, "error: prior prior far too low !!"
    assert prior < 0.50, "error: prior prior far too high !!"

    ## next, remove the priors from all scores. notabene: no prior correction on the hscore is needed.

    cor_a_corrected                    = compute_prior_away (cscore, prior)
    nscore_prior_corrected             = compute_prior_away (nscore, prior)
    fscore_prior_corrected             = compute_prior_away (fscore, prior)
    pscore_prior_corrected             = compute_prior_away (pscore, prior)
    ascore_prior_corrected             = compute_prior_away (ascore, prior)
    escore_prior_corrected             = compute_prior_away (escore, prior)
    dscore_prior_corrected             = compute_prior_away (dscore, prior)
    tscore_prior_corrected             = compute_prior_away (tscore, prior)

    ## now, do the homology correction on pscore and tscore:

    phscore_prior_corrected = pscore_prior_corrected * (1.0 - hscore)
    thscore_prior_corrected = tscore_prior_corrected * (1.0 - hscore)

    ## next, do the 1 - multiplication:

    combined_score_one_minus = (
        (1.0 - cor_a_corrected) *
        (1.0 - nscore_prior_corrected) *
        (1.0 - fscore_prior_corrected) *
        (1.0 - phscore_prior_corrected) *
        (1.0 - ascore_prior_corrected) *
        (1.0 - escore_prior_corrected) *
        (1.0 - dscore_prior_corrected) *
        (1.0 - thscore_prior_corrected)
    )

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                                 ## scale down
    combined_score += prior                                         ## and add prior.

    ## out-of-bounds check once again, to make sure:
    if combined_score < 0.0: combined_score = 0.0
    if combined_score > 1.0: combined_score = 1.0

    ## and return the result.

    return combined_score


def compute_combined_score_protein_protein_new (
        cscore,
                escore,
        dscore,
        prior = None):

    """
    >>> args = ", ".join( ["0.6"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["1.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( ["0.0"]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(0.999*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    >>> args = ", ".join( [str(1.001*local_prior_pp)]*13 )
    >>> python_perl_check("compute_combined_score_protein_protein(" + args + ")")
    """

    if prior is None:
        prior = local_prior_pp

    ## some basic checks to learn about impossible errors the hard way:

    assert prior > 0.01, "error: prior prior far too low !!"
    assert prior < 0.50, "error: prior prior far too high !!"

    ## next, remove the priors from all scores. notabene: no prior correction on the hscore is needed.

    cor_a_corrected                    = compute_prior_away (cscore, prior)
    escore_prior_corrected             = compute_prior_away (escore, prior)
    dscore_prior_corrected             = compute_prior_away (dscore, prior)


    ## next, do the 1 - multiplication:

    combined_score_one_minus = (
        (1.0 - cor_a_corrected) *
        (1.0 - escore_prior_corrected) *
        (1.0 - dscore_prior_corrected)
    )

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                                 ## scale down
    combined_score += prior                                         ## and add prior.

    ## out-of-bounds check once again, to make sure:
    if combined_score < 0.0: combined_score = 0.0
    if combined_score > 1.0: combined_score = 1.0

    ## and return the result.

    return combined_score


############################################################################################################
## subroutine: compute_combined_score_orthgroup_orthgroup ()
##
############################################################################################################

def compute_combined_score_orthgroup_orthgroup(nscore, fscore, pscore, ascore, escore, dscore, tscore, prior = None):

    if prior is None:
        prior = local_prior_pp

    ## some basic checks to learn about impossible errors the hard way:

    assert prior > 0.01, "error: prior prior far too low !!"
    assert prior < 0.50, "error: prior prior far too high !!"

    ## next, remove the priors from all scores. notabene: no prior correction on the hscore is needed.

    nscore_prior_corrected = compute_prior_away (nscore, prior)
    fscore_prior_corrected = compute_prior_away (fscore, prior)
    pscore_prior_corrected = compute_prior_away (pscore, prior)
    ascore_prior_corrected = compute_prior_away (ascore, prior)
    escore_prior_corrected = compute_prior_away (escore, prior)
    dscore_prior_corrected = compute_prior_away (dscore, prior)
    tscore_prior_corrected = compute_prior_away (tscore, prior)

    ## next, do the 1 - multiplication:

    combined_score_one_minus = (
        (1.0 - nscore_prior_corrected) *
        (1.0 - fscore_prior_corrected) *
        (1.0 - pscore_prior_corrected) *
        (1.0 - ascore_prior_corrected) *
        (1.0 - escore_prior_corrected) *
        (1.0 - dscore_prior_corrected) *
        (1.0 - tscore_prior_corrected)
    )

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*

    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                                 ## scale down
    combined_score += prior                                         ## and add prior.

    ## out-of-bounds check once again, to make sure:
    if combined_score < 0.0: combined_score = 0.0
    if combined_score > 1.0: combined_score = 1.0

    ## and return the result.

    return combined_score


############################################################################################################
## subroutine: compute_combined_score_protein_chemical ()
##
############################################################################################################

def compute_combined_score_protein_chemical(escore, escore_transferred, dscore, dscore_transferred, tscore, tscore_transferred, prior = None

    if prior is None:
        prior = local_prior_pc

    # set 0 score for channels that do not exist for protein-chemical relations
    return compute_combined_score_protein_protein(0, 0, 0, 0, 0, 0, 0, escore, escore_transferred, dscore, dscore_transferred, tscore, tscor


############################################################################################################
## subroutine: compute_combined_score_chemical_chemical ()
##
############################################################################################################

def compute_combined_score_chemical_chemical(escore, dscore, tscore, prior = None):

    if prior is None:
        prior = local_prior_cc

    # set 0 score for channels that do not exist for protein-chemical relations
    return compute_combined_score_protein_protein(0, 0, 0, 0, 0, 0, 0, escore, 0, dscore, 0, tscore, prior)




############################################################################################################
## subroutine: compute_prior_away ()
##
############################################################################################################

def compute_prior_away(score, prior):

    score = score - prior          ## substract the prior
    if score < 0: score = 0        ## avoid out-of-bounds
    score /= 1 - prior             ## correct the score so that 1.0 remains 1.0 despite prior correction
    if score > 1.0: score = 1.0    ## avoid out-of-bounds again (here: rounding/precision error).
    return score                   ## and return the corrected score.


############################################################################################################
## subroutine: combine_two_scores_generic ()
##
############################################################################################################

def combine_two_scores_generic(score1, score2, prior):

    assert prior > 0.01, "error: prior prior far too low !!"
    assert prior < 0.50, "error: prior prior far too high !!"
    assert score1 >= 0 and score2 >= 0, "error: negative scores !!"

    if score1 == 0 and score2 == 0:
        return 0

    if score1 <= prior:
        return score2

    if score2 <= prior:
        return score1

    ## first, remove the priors from both scores.
    score1_prior_corrected = compute_prior_away (score1, prior)
    score2_prior_corrected = compute_prior_away (score2, prior)

    ## next, do the 1 - multiplication:
    combined_score_one_minus = (1.0 - score1_prior_corrected) * (1.0 - score2_prior_corrected)

    ## and lastly, do the 1 - conversion again, and put back the prior *exactly once*
    combined_score = (1.0 - combined_score_one_minus)            ## 1- conversion
    combined_score *= (1.0 - prior)                                 ## scale down
    combined_score += prior                                         ## and add prior.

    ## out-of-bounds check once again, to make sure:
    if combined_score < 0.0: combined_score = 0.0
    if combined_score > 1.0: combined_score = 1.0

    ## and return the result.

    return combined_score



############################################################################################################
## subroutine: combine_two_scores_protein_protein ()
##
############################################################################################################

def combine_two_scores_protein_protein(score1, score2, prior = None):
    """
    >>> python_perl_check("combine_two_scores_protein_protein(0.5,0.5)")
    >>> python_perl_check("combine_two_scores_protein_protein(1,1)")
    >>> python_perl_check("combine_two_scores_protein_protein(0.15,0.15,0.14)")
    >>> python_perl_check("combine_two_scores_protein_protein(0.15,0.15,0.15)")
    >>> python_perl_check("combine_two_scores_protein_protein(0.15,0.15,0.151)")
    """
    if prior is None:
        prior = local_prior_pp
    return combine_two_scores_generic (score1, score2, prior)



############################################################################################################
## subroutine: combine_two_scores_protein_chemical ()
##
############################################################################################################

def combine_two_scores_protein_chemical(score1, score2, prior = None):
    """
    >>> python_perl_check("combine_two_scores_protein_chemical(0.5,0.5)")
    >>> python_perl_check("combine_two_scores_protein_chemical(1,1)")
    >>> python_perl_check("combine_two_scores_protein_chemical(0.15,0.15,0.14)")
    >>> python_perl_check("combine_two_scores_protein_chemical(0.15,0.15,0.15)")
    >>> python_perl_check("combine_two_scores_protein_chemical(0.15,0.15,0.151)")
    """
    if prior is None:
        prior = local_prior_pc
    return combine_two_scores_generic (score1, score2, prior)



############################################################################################################
## subroutine: combine_two_scores_chemical_chemical ()
##
############################################################################################################

def combine_two_scores_chemical_chemical (score1, score2, prior = None):
    """
    >>> python_perl_check("combine_two_scores_chemical_chemical(0.5,0.5)")
    >>> python_perl_check("combine_two_scores_chemical_chemical(1,1)")
    >>> python_perl_check("combine_two_scores_chemical_chemical(0.15,0.15,0.14)")
    >>> python_perl_check("combine_two_scores_chemical_chemical(0.15,0.15,0.15)")
    >>> python_perl_check("combine_two_scores_chemical_chemical(0.15,0.15,0.151)")
    """
    if prior is None:
        prior = local_prior_cc
    return combine_two_scores_generic (score1, score2, prior)


if __name__ == '__main__':
    from xmlrpclib import ServerProxy, Error

    def python_perl_check(call, server = ServerProxy("http://localhost:8047")):

        python = eval(call)
        perl = eval("server." + call)

        _python = "%.5g" % python
        _perl = "%.5g" % perl

        s =  "Call: %s\nPython result: %s\n  Perl result: %s\n" % (call, _python, _perl)
        # print >> sys.stderr, s

        assert _perl == _python, s


def main():
    import doctest
    doctest.testmod()



if __name__ == '__main__':
    main()
