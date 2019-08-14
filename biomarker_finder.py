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

class Subtype(object):
    def __init__(self, name, condition_names):
        self.name = name
        self.condition_names = condition_names
        self.conditions = []
        self.potential_biomarkers = []

    def add_data(self, df):
        self.conditions.append(df)

    def add_potential_biomarkers(self, list_of_potential_biomarkers):
        self.potential_biomarkers.append(list_of_potential_biomarkers)

class Type(object): # e.g. type of cancer. root folder
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.subtypes = self.walk_folder()

    def walk_folder(self):
        folder = self.folder_name
        subfolders = [x for x in os.walk(folder)]
        print(subfolders)
        subtypes = []
        for i in range(1, len(subfolders)):
            # name of subtype, i.e. folder with condition sheets in it 
            subtype_name = subfolders[0][1][i-1]
            print(subtype_name)
            subtype_conditions = subfolders[i][2]
            print(subtype_conditions)
            subtypes.append(Subtype(subtype_name, subtype_conditions))
        return subtypes


class BiomarkerFinder(object):
    def __init__(self, folder_name):
        self.excluded = []
        self.potential_biomarkers = []
        self.type = Type(folder_name)
        self.prepare_data()

    def prepare_data(self):
        subtypes = []
        for subtype in self.type.subtypes:
            for condition in subtype.conditions:
                spreadsheet_filename = self.type.folder_name + '/' + subtype.name + '/' + condition
                spreadsheet = self.prepare_spreadsheet(spreadsheet_filename)
                subtype.add_data(spreadsheet)
                subtypes.append(subtype)
        self.type.subtypes = subtypes

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

    def find_potential_biomarkers(self, condition_of_interest, other_conditions):
        biomarkers_in_other_conditions = [] # biomarkers to exclude are from other conditions
        for df in other_conditions:
            biomarkers_in_other_conditions.extend(list(df.index))
        potential_biomarkers = condition_of_interest[~condition_of_interest.index.isin(biomarkers_in_other_conditions)]
        shared_proteins = condition_of_interest[condition_of_interest.index.isin(biomarkers_in_other_conditions)]
        # more compilcated with list of dataframes
        for i, row in shared_proteins.iterrows():
            expr = row['up']
            # iterate over each df to compare expression, if expression is different in all, accept
            accept = True
            for df in other_conditions:
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
        for i, subtype in enumerate(self.type.subtypes):
            other_subtypes = [st for j, st in enumerate(self.type.subtypes) if j != i]
            to_compare_other_subtypes = [cond for cond in st.conditions for subtype in other_subtypes]
            for k, condition in subtype.conditions:
                if i != len(subtype.conditions) - 1:
                    to_compare = subtype.conditions[:i] + subtype.conditions[i+1:]
                else:
                    to_compare = self.data[:-1]
                to_compare.extend(to_compare_other_subtypes)
                potential_biomarkers = self.find_potential_biomarkers(condition, to_compare)
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