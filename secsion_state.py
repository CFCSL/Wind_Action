#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 14:36:14 2023

@author: namnguyen
"""

import streamlit as st
import pandas as pd
import os
from collections import defaultdict


# Initialization
if 'key' not in st.session_state:
    st.session_state['key'] = 'value'

# Session State also supports attribute based syntax
if 'key' not in st.session_state:
    st.session_state.key = 'value'

st.write(st.session_state.key)

st.session_state.key
st.session_state.key
st.session_state.key

del st.session_state['key']