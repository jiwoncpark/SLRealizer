#!/usr/bin/python2
# -*- coding: utf-8 -*-

import dropbox, sys, os
import requests

def save_file_to_dropbox(path_to_files, save_name):
    access_token = 'mwiZY0MbbvAAAAAAAAAAB3qk-8Gk9zGkaypCKGNzOeYXn0iQ1S-JHHGgee3dUru6'
    dbx = dropbox.Dropbox(access_token)
    with open(path_to_files) as f:
        dbx.files_upload(f.read(), save_name)
