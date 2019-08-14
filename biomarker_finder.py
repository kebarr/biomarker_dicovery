import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

""" Spreadsheet format: one spreadsheet per condition, i.e. early vs late, etc. conditions are labelled Group A and Group B
There can be different numbers of replicates, so need to count how many columns for group A and how many for group B
For each subtype, assuming we are doing 2 v 2 comparisons, there are 4 condition spreadsheets.
Format will be xlsx
Assume that header is always 2 lines, first line, ignore

"""

class Type(object): # e.g. type of cancer. root folder.
    def __init__(self, folder_name):
        self.folder_name = folder_name
        # may not actually need this as we're only really concerned with getting the filenames... 
        # need to trust that the folder structure and naming is correct anyway
        self.subtypes = self.walk_folder()

    def walk_folder(self):
        folder = self.folder_name
        subfolders = [x for x in os.walk(folder)]
        print(subfolders)
        subtypes = {}
        for i in range(1, len(subfolders)):
            # name of subtype, i.e. folder with condition sheets in it 
            subtype_name = subfolders[0][1][i-1]
            print(subtype_name)
            subtype_conditions = subfolders[i][2]
            print(subtype_conditions)
            subtypes[subtype_name] = subtype_conditions
        return subtypes

    def get_data_file_name(self, subtype, condition):
        subtype_name = 'Subtype ' + str(subtype)
        condition_name = subtype_name + ' Condition' + str(condition) + '.xlsx'
        return self.folder_name + '/' + subtype_name + '/' + condition_name



class BiomarkerFinder(object):
    def __init__(self, folder_name):
        self.excluded = []
        self.potential_biomarkers = []
        self.type = Type(folder_name)
        self.prepare_data()

    def prepare_data(self):
        number_subtypes = len(self.type.subtypes.keys())
        subtypes = [[] for i in range(number_subtypes)]
        for i, subtype in enumerate(self.type.subtypes.keys()):
            for spreadsheet in self.type.subtypes[subtype]:
                spreadsheet_filename = self.type.folder_name + '/' + subtype + '/' + spreadsheet
                spreadsheet = self.prepare_spreadsheet(spreadsheet_filename)
                subtypes[i].append(spreadsheet)
        self.subtypes = subtypes

    def prepare_spreadsheet(self, filename):
        # take both rows as header, then make columns manually, first 12 are from second row, afterwards, from first
        # read excel is insanely slow!!!!!
        sheet = pd.read_excel(filename, header=[1,2]) # sheet name doesn't matter as there will only be one sheet
        colnames = [sheet.columns[i][1] for i in range(12)] + [sheet.columns[i][0] + "_" + str(i) for i in range(12, len(sheet.columns))]
        sheet.columns = colnames 
        sheet = sheet.dropna() # remove columns with undefined values
        sheet = sheet[sheet['Anova (p)'] < 0.05] # filter by p value
        cond1 = sheet.filter(like="Group A", axis=1) 
        cond2 = sheet.filter(like="Group B", axis=1)
        sheet['mean1'] = cond1.mean(axis=1) # mean of cond1
        sheet['mean2'] = cond2.mean(axis=1) # mean of cond2
        filtered = sheet[~((sheet['mean1'] < 50000) & (sheet['mean2'] < 50000))]
        filtered['down'] = np.where(filtered['Highest mean condition'] == 'Group A', True, False)
        filtered['up'] = np.where(filtered['Highest mean condition'] == 'Group A', False, True)
        filtered['Accession'] = filtered['Accession'].str.split(';', n=1, expand=True)[0]
        filtered = filtered.set_index('Accession')
        return filtered

    def find_potential_biomarkers(self, df_of_potential_biomarkers, dfs_of_biomarkers_to_exclude):
        to_exclude = []
        for df in dfs_of_biomarkers_to_exclude:
            to_exclude.extend(list(df.index))
        potential_biomarkers = df_of_potential_biomarkers[~df_of_potential_biomarkers.index.isin(to_exclude)]
        shared_proteins = df_of_potential_biomarkers[df_of_potential_biomarkers.index.isin(to_exclude)]
        # more compilcated with list of dataframes
        for i, row in shared_proteins.iterrows():
            expr = row['up']
            # iterate over each df to compare expression, if expression is different in all, accept
            accept = True
            for df in dfs_of_biomarkers_to_exclude:
                try:
                    expr_other = shared_proteins.loc[i]['up']
                    if expr == expr_other:
                        # shared expression found so don't use as biomarker
                        accept = False
                        break
                except:
                    pass
            if accept == True:
                potential_biomarkers.append(row)
        return potential_biomarkers

    def find_all_potential_biomarkers(self):
        # data is all sheets
        for i, sheet in self.data:
            if i != len(self.data) - 1:
                to_compare = self.data[:i] + self.data[i+1:]
            else:
                to_compare = self.data[:-1]
            potential_biomarkers = self.find_potential_biomarkers(sheet, to_compare)
            # ones we have excluded are those indexes in sheet that are not in potential biomarkers
            # ideally want ones that are excluded based on comparison with others only, 
            excluded = sheet[~sheet.index.isin(potential_biomarkers.index)]
            self.potential_biomarkers.append(potential_biomarkers)
            self.excluded.append(excluded)

    def output(self, output_filename):
        # output potential biomarkers and excluded
        writer = pd.ExcelWriter(output_filename, engine = 'xlsxwriter')
        # for time being, just output all to one spreadsheet
        sheet_name_string = output_filename.split('.')[0]
        for i, df in enumerate(self.potential_biomarkers):
            sheet_name = 'biomarkers_' + str(i) + sheet_name_string
            df.to_excel(writer, sheet_name=sheet_name)
        for i, df in enumerate(self.excluded):
            sheet_name = 'excluded_' + str(i) + sheet_name_string
            df.to_excel(writer, sheet_name=sheet_name)
        writer.save()