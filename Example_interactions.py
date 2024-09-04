#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 01 16:01:00 2023

@author: Kabelo Mokgopa
"""
from interactions import intarna

target_seq = "CCAGAGGUUGUAACGUUGUCUAUAUAUACCCUGUAGAACCGAAUUUGUGUGGUAUCCGUAUAGUCACAGAUUCGAUUCUAGGGGAAUAUAUGGUCGAUGCAAAAACUUCA"
p= intarna('aptamers.csv', 'Aptamer',target_seq,"Pre_miR10b",True)
print(p)
