# *******************************************************************************
# Copyright (C) 2020 INAF
#
# This software is distributed under the terms of the BSD-3-Clause license
#
# Authors:
# Ambra Di Piano <ambra.dipiano@inaf.it>
# *******************************************************************************

import xml.etree.ElementTree as ET

class ManageXml():
    '''
    WRITE DOCS
    '''

    def __init__(self, xml):
        self.__xml = xml
        self.file = open(self.__xml)
        self.src_lib = ET.parse(self.file)
        self.root = self.src_lib.getroot()
        self.tsv_list = []
        self.pos = []
        self.err = []
        self.spectral = []
        self.sigma = 5
        self.default_model = True
        self.instr = 'CTA'
        self.bkg_type = 'Irf'
        self.src_att = []
        self.bkg_att = []
        self.tscalc = True
        self.if_cut = False

    # get source element ---!
    def __getSrcObj(self):
        src = self.root.findall('source')
        return src

    # skip node listed in skip element or filters ---!
    def __skipNode(self, cfg):
        src = self.__getSrcObj()
        if src.attrib[cfg.get('idAttribute')] in cfg.get('skip'):
            return True
        for filter in cfg.get('filters'):
            if src.attrib[filter.get('attribute')] == filter.get('value'):
                return True
        if len(cfg.get('selectors')) == 0:
            return False
        for select in cfg.get('selectors'):
            if src.attrib[select.get('attribute')] == select.get('value'):
                return False
        return True

    # get TS values ---!
    def getTs(self, highest=None):
        for src in self.root.findall('source'):
            if highest == None:
                if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                    tsv = float(src.attrib['ts'])
                    self.tsv_list.append(tsv)
            else:
                if src.attrib['name'] == highest:
                    tsv = float(src.attrib['ts'])
                    self.tsv_list.append(tsv)
        return self.tsv_list

    # get RA/DEC values ---!
    def getRaDec(self, highest=None):
        ra_list, dec_list = ([] for i in range(2))
        for src in self.root.findall('source'):
            if highest == None:
                if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                    ra = float(src.find('spatialModel/parameter[@name="RA"]').attrib['value'])
                    dec = float(src.find('spatialModel/parameter[@name="DEC"]').attrib['value'])
                    ra_list.append(ra)
                    dec_list.append(dec)
            else:
                if src.attrib['name'] == highest:
                    ra = float(src.find('spatialModel/parameter[@name="RA"]').attrib['value'])
                    dec = float(src.find('spatialModel/parameter[@name="DEC"]').attrib['value'])
                    ra_list.append(float(ra))
                    dec_list.append(float(dec))
        self.pos = [ra_list, dec_list]
        return self.pos

    # get Gaussian sigma values ---!
    def getConfInt(self, highest=None):
        ra_list, dec_list = ([] for i in range(2))
        for src in self.root.findall('source'):
            if highest == None:
                if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                    ra = float(src.find('spatialModel/parameter[@name="RA"]').attrib['value'])
                    dec = float(src.find('spatialModel/parameter[@name="DEC"]').attrib['value'])
                    ra_list.append(ra)
                    dec_list.append(dec)
            else:
                if src.attrib['name'] == highest:
                    ra = float(src.find('spatialModel/parameter[@name="RA"]').attrib['value'])
                    dec = float(src.find('spatialModel/parameter[@name="DEC"]').attrib['value'])
                    ra_list.append(ra)
                    dec_list.append(dec)
        self.err = [ra_list, dec_list]
        return self.err

    # get spectral parameter values ---1
    def getSpectral(self, highest=None):
        index_list, pref_list, pivot_list = ([] for i in range(3))
        if self.if_cut is True :
            cutoff_list = []

        for src in self.root.findall('source'):
            if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                if highest == None:
                    index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
                        src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
                    pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
                        src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
                    pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
                        src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
                    index_list.append(index)
                    pref_list.append(pref)
                    pivot_list.append(pivot)
                    if self.if_cut is True :
                        cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                            src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
                        cutoff_list.append(cutoff)
                else:
                    if src.attrib['name'] == highest:
                        index = float(src.find('spectrum/parameter[@name="Index"]').attrib['value']) * float(
                            src.find('spectrum/parameter[@name="Index"]').attrib['scale'])
                        pref = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['value']) * float(
                            src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
                        pivot = float(src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['value']) * float(
                            src.find('spectrum/parameter[@name="PivotEnergy"]').attrib['scale'])
                        index_list.append(index)
                        pref_list.append(pref)
                        pivot_list.append(pivot)
                        if self.if_cut is True:
                            cutoff = float(src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['value']) * float(
                                src.find('spectrum/parameter[@name="CutoffEnergy"]').attrib['scale'])
                            cutoff_list.append(cutoff)

        if self.if_cut is False:
            self.spectral = [index_list, pref_list, pivot_list]
        else:
            self.spectral = [index_list, pref_list, pivot_list, cutoff_list]
        return self.spectral

    # get prefact error values ---!
    def getPrefError(self, highest=None):
        err_list = []
        for src in self.root.findall('source'):
            if highest == None:
                if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                    err = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['error']) * float(
                        src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
                    err_list.append(err)
            else:
                if src.attrib['name'] == highest:
                    err = float(src.find('spectrum/parameter[@name="Prefactor"]').attrib['error']) * float(
                        src.find('spectrum/parameter[@name="Prefactor"]').attrib['scale'])
                    err_list.append(err)
        self.err = err_list
        return self.err

    # save xml files ---!
    def __saveXml(self):
        self.src_lib.write(self.__xml, encoding="UTF-8", xml_declaration=True)
        return

    # set a default spectral model and bkg ---!
    def __setModel(self):
        if self.default_model is True:
            att_prefactor = {'name': 'Prefactor', 'scale': '1e-16', 'value': '5.7', 'min': '1e-07', 'max': '1e7', 'free': '1'}
            att_index = {'name': 'Index', 'scale': '-1', 'value': '2.48', 'min': '0', 'max': '5.0', 'free': '1'}
            att_pivot = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1.0', 'min': '1e-07', 'max': '1000.0', 'free': '0'}
            bkg_prefactor = {'name': 'Prefactor', 'scale': '1', 'value': '1.0', 'min': '1e-03', 'max': '1e+3.0', 'free': '1'}
            bkg_index = {'name': 'Index', 'scale': '1.0', 'value': '0.0', 'min': '-5', 'max': '+5.0', 'free': '1'}
            bkg_pivot = {'name': 'PivotEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '0'}

            self.src_att = [att_prefactor, att_index, att_pivot]
            self.bkg_att = [bkg_prefactor, bkg_index, bkg_pivot]
            if self.if_cut is True:
                att_cutoff = {'name': 'CutoffEnergy', 'scale': '1e6', 'value': '1.0', 'min': '0.01', 'max': '1000.0', 'free': '1'}
                self.src_att.append(att_cutoff)
        else:
            pass

    # set tscalc to 1 ---!
    def setTsTrue(self):
        for src in self.root.findall('source'):
            if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                src.set('tscalc', '1')
        self.__saveXml()
        return

    # modeify the spectral component of candidate list ---!
    def modXml(self, overwrite=True):
        self.__setModel()
        # source ---!
        i = 0
        for src in self.root.findall('source'):
            i += 1
            if src.attrib['name'] not in ('Background', 'CTABackgroundModel'):
                src.set('tscalc', '1') if self.tscalc is True else None
                # remove spectral component ---!
                rm = src.find('spectrum')
                src.remove(rm)
                # new spectrum ---!
                if self.if_cut:
                    spc = ET.SubElement(src, 'spectrum', attrib={'type': 'ExponentialCutoffPowerLaw'})
                else:
                    spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
                spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
                spc.tail = '\n\t\t'.replace('\t', ' ' * 2)
                #src.insert(0, spc)
                # new spectral params ---!
                for j in range(len(self.src_att)):
                    prm = ET.SubElement(spc, 'parameter', attrib=self.src_att[j])
                    if prm.attrib['name'] == 'Prefactor' and i > 1:
                        prm.set('value', str(float(prm.attrib['value']) / 2 ** (i - 1)))
                    prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.src_att) else '\n\t\t'.replace('\t', ' ' * 2)
                    #spc.insert(j, prm)
            # background ---!
            else:
                # set bkg attributes ---!
                src.set('instrument', '%s' % self.instr.upper()) if self.instr.capitalize() != 'None' else None
                if self.bkg_type.capitalize() == 'Aeff' or self.bkg_type.capitalize() == 'Irf':
                    src.set('type', 'CTA%sBackground' % self.bkg_type.capitalize())
                if self.bkg_type.capitalize() == 'Racc':
                    src.set('type', 'RadialAcceptance')
                # remove spectral component ---!
                rm = src.find('spectrum')
                src.remove(rm)
                # new bkg spectrum ---!
                spc = ET.SubElement(src, 'spectrum', attrib={'type': 'PowerLaw'})
                spc.text = '\n\t\t\t'.replace('\t', ' ' * 2)
                spc.tail = '\n\t'.replace('\t', ' ' * 2)
                # new bkg params ---!
                for j in range(len(self.bkg_att)):
                    prm = ET.SubElement(spc, 'parameter', attrib=self.bkg_att[j])
                    prm.tail = '\n\t\t\t'.replace('\t', ' ' * 2) if j < len(self.bkg_att) else '\n\t\t'.replace('\t', ' ' * 2)

        # instead of override original xml, save to a new one with suffix "_mod" ---!
        if not overwrite:
            self.__xml = self.__xml.replace('.xml', '_mod.xml')
        self.__saveXml()
        return

    # free and fix parameters for max like computation ---!
    def parametersFreeFixed(self, src_free=['Prefactor'], bkg_free=['Prefactor', 'Index']):
        for src in self.root.findall('source'):
            if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                for prm in src.findall('*/parameter'):
                    if prm.attrib['name'] not in src_free:
                        prm.set('free', '0')
                    else:
                        prm.set('free', '1')
            else:
                for prm in src.findall('*/parameter'):
                    if prm.attrib['name'] not in bkg_free:
                        prm.set('free', '0')
                    else:
                        prm.set('free', '1')

        self.setTsTrue() if self.tscalc is True else None
        self.__saveXml()
        return

    # sort candidates by their ts value ---!
    def sortSrcTs(self):
        src = self.root.findall("*[@ts]")
        self.root[:-1] = sorted(src, key=lambda el: (el.tag, el.attrib['ts']), reverse=True)
        from_highest = []
        for src in self.root.findall("*[@ts]"):
            from_highest.append(src.attrib['name'])
        self.__saveXml()
        if len(from_highest) == 0:
            from_highest = [None]
        return from_highest

    # close xml ---!
    def closeXml(self):
        self.file.close()
        return

    # find observation filename ---!
    def getRunList(self):
        filenames = []
        for obs in self.root.findall('observation/parameter[@name="EventList"]'):
            filenames.append(obs.attrib['file'])
        return filenames

    def setModelParameters(self, source, parameters=(), values=()):
        parameters = tuple(parameters)
        values = tuple(values)
        for src in self.root.findall('source'):
            if src.attrib['name'] != 'Background' and src.attrib['name'] != 'CTABackgroundModel':
                src.set('name', source)
                for i, prm in enumerate(parameters):
                    src.find('*/parameter[@name="%s"]' % prm).set('value', str(values[i]))
        self.__saveXml()
        return