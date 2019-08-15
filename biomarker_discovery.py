import click # command liine options

from biomarker_finder import BiomarkerFinder
# source ~/python3/bin/activate - for python 3
# for now write as command line with input and output names specified

# in input option, nargs is number of input files, we have 4. type tells it its a filepath.
# actually just take folder name 
@click.command() 
@click.option("--input_folder", "-i", help="Folder containing input data, expected to contain folders called 'Subtype1' etc, which then contain spreadsheets with names 'Subtype1 Condition1' etc", required=True) 
@click.option("--output", "-o", help="Output filename, if you want a different name/location than the default")
@click.option("--condition1", "-c", help="Condition you want to find biomarkers for", default="Condition1")
@click.option("--condition2", "-k", help="Other condition in subtype that you want to compare against", default="Condition2")
@click.option("--subtype", "-s", help="Subtype you want to find biomarkers for", default="Subtype1")
@click.option("--flowchart", "-f", help="flowchart logic, can either be basic (default), comparing between two conditions in one subtype, or 'all', which additionally compares the other subtypes", type=click.Choice(['basic', 'all']), default="basic")
def run_analysis(input_folder, condition1, condition2, subtype, flowchart, output=None):
    """Run the biomarker discovery analysis"""
    print("intitalising data...")
    bf = BiomarkerFinder(input_folder)
    print("data initialised, running flowchart %s" % (flowchart))
    if flowchart == "basic":
        bf.compare_two_conditions_in_same_subtype(subtype_name=subtype, condition_name1=condition1, condition_name2=condition2, out_filename=output)
    elif flowchart == "all":
        bf.find_diagnosis_biomarkers(subtype_name=subtype, condition_name1=condition1, condition_name2=condition2, out_filename=output)

# don't think adding button to excel will work so try anf make an executable with a little GUI
# https://docs.python-guide.org/shipping/freezing/
if __name__ == '__main__':
    run_analysis()