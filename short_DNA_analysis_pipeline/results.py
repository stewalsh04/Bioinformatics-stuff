# -*- coding: utf-8 -*-
"""
Class to store the results from the pipeline
"""

class ngsRseults:
    def __init__(self, b64matches:str, b64mismatches:str, b64totalreads:str):
        self.b64mtches = b64matches # base64 match figure
        self.b64mismatches = b64mismatches # base64 mis-matches figure
        self.b64totalreads = b64totalreads # base64 total reads figure 
