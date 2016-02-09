#!/usr/bin/env python
"""Validate the arguments that are passed. The functions defined here will check
that files are readable/writeable, email addresses are valid, and that user
specifed values are in the expected ranges."""

import os
import re
from Bio import SeqIO


def valid_dir(dname):
    """Return True if a given directory exists and has read/write permissions,
    return False otherwise."""
    #   First, expand any tildes that may be in the path
    fpath = os.path.expanduser(dname)
    #   Then resolve any relative paths that might have been given
    abspath = os.path.abspath(fpath)
    #   os.path.isdir() checks whether a given path is a directory or not
    #   (as opposed to a file or a device)
    if os.path.isdir(abspath):
        try:
            ls = os.listdir(abspath)
            return True
        except (OSError, PermissionError, FileNotFoundError):
            return False
    else:
        return False


def valid_email(email):
    """Returns True if a specified email address is valid, and False if not.
    This doesn't do a perfect check, as the official specification is very
    complicated, but this simple regex will do for now."""
    #   The regex is as follows:
    #       1 or more non-whitespace characters
    #       @
    #       1 or more of ([1 or more non-whitespace characters].)
    #       1 non-whitespace character
    email_regex = '\S+@(\S+\.)+\S'
    if re.match(email_regex, email):
        return True
    else:
        return False


def valid_fasta(fasta):
    """Returns True if a specified FASTA file is valid and False if not. Catches
    an exception thrown by parsing the provided file and converting to a list of
    sequences."""
    try:
        seq_list = list(SeqIO.parse(fasta, 'fasta'))
        for s in seq_list:
            #   If any sequence of the FASTA file has illegal characters, we
            #   raise and error so that it will get caught as an invalid file
            if re.match('[^ATCGMRWSKYBDHVN\[\/\]]', str(s.seq)):
                raise ValueError
        return True
    except:
        #   This will catch cases where the file provided is totally different
        #   from a FASTA file, or contains weird characters.
        return False


def validate_args(args):
    """Validate the arguments. Files should be readable, directories should have
    read access, FASTA files should be valid, and if supplied, the contextual
    sequence length should be a positive integer. Return codes are error codes
    and mean the following:
        NoneType) No errors.
        1) Required arguments not supplied. Print usage.
        2) Email is not valid.
        3) Context length is not valid.
        4) Directory is not valid.
        5) FASTA is not valid."""
    #   Convert the arguments to a dictionary
    arg_dict = args.__dict__
    #   Check our required argument groups:
    #       1) Fasta/Directory
    #       2) Context length/GBS/Illumina
    #       3) Email
    if not (arg_dict['fasta_file'] or arg_dict['dir']):
        return 1
    elif not (arg_dict['clen'] or arg_dict['gbs'] or arg_dict['illumina']):
        return 1
    elif not arg_dict['email']:
        return 1

    #   Then check the email
    if not valid_email(arg_dict['email']):
        return 2

    #   If the context length is not valid, return a flag that says to fix the
    #   length
    if arg_dict['clen']:
        try:
            c = int(arg_dict['clen'])
            assert c > 0, 'Negative'
        except (TypeError, AssertionError):
            return 3

    #   If the directory is not valid, do the same as above
    if arg_dict['dir']:
        if not valid_dir(arg_dict['dir']):
            return 4

    #   Then, the FASTA file
    elif arg_dict['fasta_file']:
        if not valid_fasta(arg_dict['fasta_file']):
            return 5
    return

