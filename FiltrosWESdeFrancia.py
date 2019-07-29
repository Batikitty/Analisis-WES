# -*- coding: utf-8 -*-
"""
Author: Manuela Molina Ospina
Date: 27 Julio 2019
Utilidad: Filtrar las tablas de WES originadas a partir del pipeline de Francia
"""
import pandas as pd

#Leer csv
dt = pd.read_csv('/Volumes/lymphocytex/Documents/CSV JL 20190709/COL1533.V4.2.csv', header=0)

#Filtrar por frecuencias menores o iguales al 100%
dt['gnomADGenomes_AF'] = pd.to_numeric(dt['gnomADGenomes_AF'], errors = 'coerce')
dt['gnomADGenomes_AF'] = dt['gnomADGenomes_AF'].fillna(0)
freqMen = dt['gnomADGenomes_AF'] <= 0.01
#freqMay = dt['gnomADGenomes_AF'] > 0.01

#Filtrar por Zigocidad
fHom = dt['zygo'] == 'hom'
fHet = dt['zygo'] == 'het'
#fNoHom = dt['zygo'] != 'hom'
#fNoHet = dt['zygo'] != 'het'

#Quitar regiones no codificantes
fFun = (dt['function'] != 'intron') & (dt['function'] != 'nc_mirna') & (dt['function'] != 'upstream') & (dt['function'] != 'utr_3') & (dt['function'] != 'utr_5') & (dt['function'] != 'downstream')
#& (dt['function'] != 'synonymous') Se están dejando las sinónimas
#fNoFun = (dt['function'] == 'intron') | (dt['function'] == 'nc_mirna')

#Filtrar por GDI, ideal para quitar todo lo mayor a 15, entre el número es más alto más mutaciones hay en ese gen en gente sana
dt['GDI'] = pd.to_numeric(dt['GDI'], errors = 'coerce')
dt['GDI'] = dt['GDI'].fillna(0)
fGDIMen = dt['GDI'] <= 15.0
#fGDIMay = dt['GDI'] > 15.0

#Quitar el MSC low porque no predice alto daño
#fMSC = dt['MSC_99%_Pred'] != 'HIGH'
fMSCHigh = dt['MSC_99%_Pred'] != 'LOW' # Toco corregir porque en algunos habían espacios en blanco

#¿está en la black list o no?
#fBLy = dt['BL'] == 'yes'
fBLn = dt['BL'] == 'no'

#Filtrar por número de reads que cubren dicha zona
dt['DP'] = pd.to_numeric(dt['DP'], errors = 'coerce')
dt['DP'] = dt['DP'].fillna(0)
#fDPMen = dt['DP'] <= 2
fDPMay = dt['DP'] > 2

#¿Está en la lista de PID?
fPID = dt['KnownPID'] == 'YES'
fNoPID = dt['KnownPID'] != 'YES'

#¿Está en la lista de predict PID?
fPPID = dt['PredictedPID'] == 'YES'
fNoPPID = dt['PredictedPID'] != 'YES'


# Datos para plantilla de Carlos hacer homo y het por separado para evitar confusiones
'''
#Homos
d2 = dt[freqMen]
d3 = dt[freqMen & fHom]
d4 = dt[freqMen & fHom & fFun]
d5 = dt[freqMen & fHom & fFun & fGDIMen]
d6 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh]
d7 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn]
d8 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay]
t1 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fPID]
t2 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID]
t3 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID & fPPID]
t4 = dt[freqMen & fHom & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID & fNoPPID]
'''
#Heteros
d2 = dt[freqMen]
d3 = dt[freqMen & fHet]
d4 = dt[freqMen & fHet & fFun]
d5 = dt[freqMen & fHet & fFun & fGDIMen]
d6 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh]
d7 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn]
d8 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay]
t1 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fPID]
t2 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID]
t3 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID & fPPID]
t4 = dt[freqMen & fHet & fFun & fGDIMen  & fMSCHigh & fBLn & fDPMay & fNoPID & fNoPPID]

#Cambiar nombre del archivo al antojo y ubicación
export_csv = t1.to_csv (r'/Volumes/lymphocytex/Documents/Archivos LUCIA /Analisis IDP/HetPID.csv', index = None, header=True)
export_csv = t3.to_csv (r'/Volumes/lymphocytex/Documents/Archivos LUCIA /Analisis IDP/HetPredictPID.csv', index = None, header=True)
export_csv = t4.to_csv (r'/Volumes/lymphocytex/Documents/Archivos LUCIA /Analisis IDP/HetNoPredictPID.csv', index = None, header=True)

def Carlos (dates):
    for x in dates:
        print (x['chr'].count())

filters = [dt, d2, d3, d4, d5, d6, d7, d8, t1, t2, t3, t4]
Carlos (filters)
