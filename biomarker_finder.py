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
        self.potential_biomarkers = {}

    def add_data(self, df):
        self.conditions.append(df)

    def add_potential_biomarkers(self, condition_name, list_of_potential_biomarkers):
        self.potential_biomarkers[condition_name] = list_of_potential_biomarkers

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
        subtypes = []
        print(subfolders)
        for i in range(len(subfolders[0][1])):
            # name of subtype, i.e. folder with condition sheets in it 
            subtype_name = subfolders[0][1][i]
            for j in range(1, len(subfolders)):
                if subtype_name in subfolders[j][0]:
                    print("subtype: ", str(subtype_name))
                    if not os.path.exists(folder + '/' + subtype_name + '/results'):
                        print("making dir ", folder + '/' + subtype_name  + '/results')
                        os.mkdir(folder + '/' + subtype_name + '/results')
                    subtype_conditions = [x.split(' ')[1].split('.')[0] for x in subfolders[j][2]]
                    print('conditions: ', subtype_conditions)
                    subtypes.append(Subtype(subtype_name, subtype_conditions))
                    break
        return subtypes

    def get_subtype(self, subtype_name):
        for subtype in self.subtypes:
            if subtype.name == subtype_name:
                return subtype
        raise ValueError("No subtype called %s" % (subtype_name))

class BiomarkerFinder(object):
    def __init__(self, folder_name):
        self.discarded = {}
        self.potential_biomarkers = []
        self.type = Type(folder_name)
        self.prepare_data()

    def prepare_data(self):
        for subtype in self.type.subtypes:
            print("preparing %s data" % subtype.name)
            for condition in subtype.condition_names:
                spreadsheet_filename = self.type.folder_name + '/' + subtype.name + '/' + subtype.name + ' ' + condition + '.xlsx'
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
        print("spreadsheet %s contains %d entries after cutoffs" % (filename, len(filtered)))
        filtered['up/down'] = np.where(filtered['Highest mean condition'] == 'Group A', 'down', 'up')
        filtered['Accession'] = filtered['Accession'].str.split(';', n=1, expand=True)[0]
        filtered = filtered.set_index('Accession')
        return filtered

    
    def find_potential_biomarkers(self, condition_of_interest, other_conditions, discarded=pd.DataFrame()):
        biomarkers_in_other_conditions = [] 
        for df in other_conditions:
            biomarkers_in_other_conditions.extend(list(df.index))
        potential_biomarkers = condition_of_interest[~condition_of_interest.index.isin(biomarkers_in_other_conditions)]
        print("found %d proteins not present in other condition" % (len(potential_biomarkers)))
        shared_proteins = condition_of_interest[condition_of_interest.index.isin(biomarkers_in_other_conditions)]
        print("%d proteins shared between condition of interest and others" % (len(shared_proteins)))
        for i, row in shared_proteins.iterrows():
            expr = row['up/down']
            # iterate over each df to compare expression, if expression is different in all, accept
            accept = True
            for df in other_conditions:
                try:
                    expr_other = df.loc[i]['up/down']
                    if expr == expr_other:
                        # shared expression found so don't use as biomarker
                        accept = False
                        discarded = discarded.append(condition_of_interest.loc[i], sort=False)
                        break
                except:
                    pass
            if accept == True:
                potential_biomarkers = potential_biomarkers.append(row, sort=False)
        return potential_biomarkers, discarded

    def compare_two_conditions_in_same_subtype(self, subtype_name='Subtype1', condition_name1='Condition1', condition_name2='Condition2', only=False, out_filename=None):
        subtype = self.type.get_subtype(subtype_name)
        self.discarded[subtype_name + condition_name1] = {}
        for i in range(len(subtype.conditions)):
            if subtype.condition_names[i] == condition_name1:
                condition = subtype.conditions[i]
            elif subtype.condition_names[i] == condition_name2:
                other_condition = subtype.conditions[i]
        potential_biomarkers, discarded = self.find_potential_biomarkers(condition, [other_condition])
        self.discarded[subtype_name + condition_name1][condition_name2] = discarded 
        if only:
            if out_filename:
                self.write_potential_biomarkers_to_file(subtype_name, condition_name1, potential_biomarkers, out_filename)
            else:
                self.write_potential_biomarkers_to_file(subtype_name, condition_name1, potential_biomarkers)
            print("found %d potential biomarkers " % (len(potential_biomarkers)))
        return potential_biomarkers

    def find_diagnosis_biomarkers(self, subtype_name='Subtype1', condition_name1='Condition1', condition_name2='Condition2', other_subtypes=['Subtype2', 'Subtype3'], other_conditions=['Condition1', 'Condition2', 'Condition3'], out_filename=None):
        if subtype_name in other_subtypes:
            raise ValueError("Subtype to test %s in subtypes to compare, cannot compare against itself" % subtype_name)
        subtype = self.type.get_subtype(subtype_name)
        # from flow diagram- Control is group A - which corresponds to meanA 
        # potential biomarkers = condition vs control; is it in condition 2 vs control? if so, is expression the same?
        # if passed for that subtype, compare to conditions 1 to 3 of other subtypes in the same way. 
        potential_biomarkers = self.compare_two_conditions_in_same_subtype(subtype_name=subtype_name, condition_name1=condition_name1, condition_name2=condition_name2, out_filename=out_filename)
        to_compare = []
        # now check these against subtypes 2 and 3, conditions 1-3
        for st_name in other_subtypes:
            st = self.type.get_subtype(st_name)
            print("comparing %s to %s" % (subtype_name, st_name))
            for i, condition in enumerate(st.conditions):
                if st.condition_names[i] in other_conditions:
                    to_compare.append(condition)
        potential_biomarkers, discarded = self.find_potential_biomarkers(potential_biomarkers, to_compare) 
        print("found %d potential diagnosis biomarkers " % (len(potential_biomarkers)))
        subtype.add_potential_biomarkers(condition_name1, potential_biomarkers)
        self.write_potential_biomarkers_to_file(subtype_name, condition_name1, potential_biomarkers, out_filename)
        discarded_names = list(discarded.index)
        print(len(discarded_names), discarded_names)
        self.output_discarded(other_subtypes, other_conditions, out_filename)
        return potential_biomarkers


    def find_potential_monitoring_biomarkers(self, subtype_name='Subtype1', condition_name='Condition4', other_subtypes=['Subtype2', 'Subtype3'], other_conditions=['Condition1', 'Condition2', 'Condition3'], out_filename=None):
        if subtype_name in other_subtypes:
            raise ValueError("Subtype to test %s in subtypes to compare, cannot compare against itself" % subtype_name)
        subtype = self.type.get_subtype(subtype_name)
        for i in range(len(subtype.conditions)):
            if subtype.condition_names[i] == condition_name:
                condition = subtype.conditions[i]
        to_compare = []
        # now check these against subtypes 2 and 3, conditions 1-3
        for st_name in other_subtypes:
            st = self.type.get_subtype(st_name)
            print("comparing %s to %s" % (subtype_name, st_name))
            for i, c in enumerate(st.conditions):
                if st.condition_names[i] in other_conditions:
                    to_compare.append(c)
        potential_biomarkers, _ = self.find_potential_biomarkers(condition, to_compare) 
        subtype.add_potential_biomarkers(condition_name, potential_biomarkers)
        if len(potential_biomarkers) > 0:
            self.write_potential_biomarkers_to_file(subtype_name, condition_name, potential_biomarkers, out_filename)

    def output_discarded(self, other_subtypes, other_conditions, out_filename=None):
        if out_filename:
            discarded_filename = 'discarded_' + out_filename
        else:
            discarded_filename = 'discarded.txt'
        with open(discarded_filename, 'a+') as f:
            for key, value in self.discarded.items():
                f.write(key + "\n") # this is discarded based on comparison between key (subtype + cond) and condition within same subtype (key2)
                # need to exclude ones that would be discarded based on another subtype, so need to search other subtypes and conditions and see if it would also have been excluded based on that, and if so, do not output
                for key2, value2 in value.items():
                    print("################################")
                    print(len(value2.index), value2.index)
                    f.write(key2 + "\n")
                    final = self.keep_in_discarded_output(value2, other_subtypes, other_conditions)
                    # TODO: these are values in one we're not interested in, want them from base comparison
                    final['Gene name'] = final['Description'].str.split('GN=', expand=True)[1].str.split(" PE=", expand=True)[0] # i will likely go to hell for this
                    final['Log2 fold change'] = np.log2(final['meanB']) - np.log2(final['meanA'])
                    final['Log10 fold change'] = np.log10(final['meanB']) - np.log10(final['meanA'])
                    final.to_csv(f)

    def keep_in_discarded_output(self, discarded, other_subtypes, other_conditions):
        to_check = pd.DataFrame()
        for t in other_subtypes:
            st = self.type.get_subtype(t)
            for c in other_conditions:
                c_st = st.get_condition(c)
                to_check =  to_check.append(c_st, sort=False)
        # need inner join on index, then ones with equal highest mean condition
        in_common = discarded.join(to_check, how='inner', lsuffix='_A', rsuffix='_B')
        print(len(in_common.index))
        to_exclude = in_common[in_common['Highest mean condition_A'] == in_common['Highest mean condition_B']].index
        print(len(to_exclude))
        # to exclude should be ones that are also in other conditions, 
        final = discarded[~discarded.index.isin(to_exclude)]
        return final

    def write_potential_biomarkers_to_file(self, subtype_name, condition_name, potential_biomarkers, out_filename=None):
        if not out_filename:
            out_filename = self.type.folder_name + '/' + subtype_name + '/results/' + subtype_name + ' ' + condition_name + '_potential_biomarkers.csv'
        potential_biomarkers['Gene name'] = potential_biomarkers['Description'].str.split('GN=', expand=True)[1].str.split(" PE=", expand=True)[0] # i will likely go to hell for this
        potential_biomarkers['Description'] = potential_biomarkers['Description'].str.split('OS', 0, expand=True)[0] # get everything before 'OS'
        potential_biomarkers['Log2 fold change'] = np.log2(potential_biomarkers['meanB']) - np.log2(potential_biomarkers['meanA'])
        potential_biomarkers['Log10 fold change'] = np.log10(potential_biomarkers['meanB']) - np.log10(potential_biomarkers['meanA'])
        potential_biomarkers['Anova (p)'].replace({0:0.000000001}) # hack to avoid infinite values
        potential_biomarkers['-log 10 p'] = -(np.log10(potential_biomarkers['Anova (p)']))
        out_df = potential_biomarkers[['Gene name', 'Max fold change', 'Log2 fold change', 'Anova (p)', '-log 10 p', 'Description', 'Highest mean condition']]
        out_df = out_df.set_index('Gene name', drop=False)
        out_df.to_csv(out_filename)

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

    # https://github.com/reneshbedre/bioinfokit/blob/master/bioinfokit/visuz.py
    # actually just steal his function so i can add this to class without saving data to file and reopening
    def volcano(self, d, lfc="Log2 fold change", pv="Anova (p)", lfc_thr=2, pv_thr=0.05, plot_name='volcano.png', plot_title = "Subtype1 Condition1"):
        # load csv data file
        d.loc[(d[lfc] >= 2) & (d[pv] < pv_thr), 'color'] = "green"  # upregulated
        d.loc[(d[lfc] <= -2) & (d[pv] < pv_thr), 'color'] = "red"  # downregulated
        d['color'].fillna('grey', inplace=True)  # intermediate
        d[pv].replace({0:0.000000001}) # hack to avoid infinite values
        d['logpv'] = -(np.log10(d[pv]))
        # plot
        plt.scatter(d[lfc], d['logpv'], c=d['color'])
        plt.xlabel('log2 Fold Change', fontsize=15, fontname="sans-serif", fontweight="bold")
        plt.ylabel('-log10(P-value)', fontsize=15, fontname="sans-serif", fontweight="bold")
        plt.xticks(fontsize=12, fontname="sans-serif")
        plt.yticks(fontsize=12, fontname="sans-serif")
        plt.title(plot_title)
        plt.savefig(plot_name, format='png', bbox_inches='tight', dpi=300)
        plt.close()


# need ones that would have been accepted based on comparisons with other subtypes, but are discarded due to comparison with condition 2
# so make list of those discarded in other subtypes, and take those from condition2 that don't appear in this list
# its already been discarded so, not in shared indices- so need to look in original data for each subtype.
# so for each discarded row, need to find out if that row is present in other subtypes and conditions tested, and if so, whether the expression is the same
