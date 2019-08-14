import click # command liine options

from biomarker_finder import BiomarkerFinder
# source ~/python3/bin/activate - for python 3
# for now write as command line with input and output names specified

# in input option, nargs is number of input files, we have 4. type tells it its a filepath.
# actually just take folder name 
@click.command() 
@click.option("--input folder", "-i", help="Folder containing input data, expected to contain folders called 'Subtype1' etc, which then contain spreadsheets with names 'Subtype1 Condition1' etc", required=True) 
@click.option("--output", "-o", help="Output filename, if you want a different name/location than the default")
@click.option("--condition1", "-c", help="Condition you want to find biomarkers for", default="Condition1")
@click.option("--condition2", "-k", help="Other condition in subtype that you want to compare against", default="Condition2")
@click.option("--subtype", "-subtype", help="Subtype you want to find biomarkers for", default="Subtype1")
def run_analysis(input_folder, condition1, condition2, subtype, output_filename=None):
    """Run the biomarker discovery analysis"""
    bf = BiomarkerFinder(input_folder)
    bf.find_diagnosis_biomarkers(subtype_name=subtype, condition_name1=condition1, condition_name2=condition2, out_filename=output_filename)


if __name__ == '__main__':
    run_analysis()