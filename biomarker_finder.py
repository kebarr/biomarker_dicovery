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
        """Holder for subtype data. 
            Assumes that ordering of condition names, conditions, and potential biomarkers are always the same.
        """
        self.name = name
        self.condition_names = condition_names
        self.conditions = []
        self.potential_biomarkers = []

    def add_data(self, df):
        self.conditions.append(df)

    def add_potential_biomarkers(self, list_of_potential_biomarkers):
        self.potential_biomarkers.append(list_of_potential_biomarkers)

    def get_condition(self, condition_name):
        index = self.condition_names.index(condition_name)
        return self.conditions[index]

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

    def get_subtype(self, subtype_name):
        for subtype in self.subtypes:
            if subtype.name == subtype_name:
                return subtype
        raise ValueError("No subtype called %s" % (subtype_name))

class BiomarkerFinder(object):
    def __init__(self, folder_name):
        self.excluded = []
        self.potential_biomarkers = []
        self.type = Type(folder_name)
        self.prepare_data()

    def prepare_data(self):
        for subtype in self.type.subtypes:
            print(subtype.name)
            for condition in subtype.condition_names:
                spreadsheet_filename = self.type.folder_name + '/' + subtype.name + '/' + condition
                spreadsheet = self.prepare_spreadsheet(spreadsheet_filename)
                subtype.add_data(spreadsheet)

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
        sheet['meanA'] = cond1.mean(axis=1) # mean of cond1
        sheet['meanB'] = cond2.mean(axis=1) # mean of cond2
        filtered = sheet[~((sheet['meanA'] < 50000) & (sheet['meanB'] < 50000))]
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
        print("len potential biomarkers: %d" % (len(potential_biomarkers)))
        shared_proteins = condition_of_interest[condition_of_interest.index.isin(biomarkers_in_other_conditions)]
        print("len shared_proteins: %d" % (len(shared_proteins)))
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

    def find_diagnosis_biomarkers(self, subtype_name='Subtype1', condition_name1='Condition1', condition_name2='Condition2'):
        # from flow diagram- Control is group A - which corresponds to meanA 
        # potential biomarkers = condition vs control; is it in condition 2 vs control? if so, is expression the same?
        # if passed for that subtype, compare to conditions 1 to 3 of other subtypes in the same way. 
        subtype = self.type.get_subtype(subtype_name)
        condition_sheet_name1 = subtype_name + ' ' + condition_name1 + '.xlsx'
        condition_sheet_name2 = subtype_name + ' ' + condition_name2 + '.xlsx'
        for i in range(len(subtype.conditions)):
            print(subtype.condition_names[i])
            if subtype.condition_names[i] == condition_sheet_name1:
                condition = subtype.conditions[i]
            elif subtype.condition_names[i] == condition_sheet_name2:
                other_condition = subtype.conditions[i]
        potential_biomarkers = self.find_potential_biomarkers(condition, [other_condition])
        print(len(potential_biomarkers))


    # this doesn't actually represent flow diagram... it just does all against all comparison
    def find_all_potential_biomarkers(self):
        for i, subtype in enumerate(self.type.subtypes):
            other_subtypes = [st for j, st in enumerate(self.type.subtypes) if j != i]
            print("len other subtypes: %d " % len(other_subtypes))
            to_compare_other_subtypes = [cond for st in other_subtypes for cond in st.conditions]
            for k, condition in enumerate(subtype.conditions):
                print('k ', k, ' len condition: ', len(condition))
                if k != len(subtype.conditions) - 1:
                    to_compare = subtype.conditions[:k] + subtype.conditions[k+1:]
                else:
                    to_compare = subtype.conditions[:-1]
                to_compare.extend(to_compare_other_subtypes)
                print(len(to_compare))
                potential_biomarkers = self.find_potential_biomarkers(condition, to_compare)
                # ones we have excluded are those indexes in sheet that are not in potential biomarkers
                # ideally want ones that are excluded based on comparison with others only, 
                excluded = condition[~condition.index.isin(potential_biomarkers.index)]
                subtype.potential_biomarkers.append(potential_biomarkers)
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