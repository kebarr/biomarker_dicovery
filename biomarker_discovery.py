import click # command liine options

from biomarker_finder import BiomarkerFinder
# source ~/python3/bin/activate - for python 3
# for now write as command line with input and output names specified

# in input option, nargs is number of input files, we have 4. type tells it its a filepath.
# actually just take folder name 
@click.command() 
@click.option("--input folder", "-i", help="Folder containing input data, expected to contain folders called 'Subtype 1' etc, which then contain spreadsheets with naes 'Subtype 1 Condition 1' etc", nargs=4, type=click.Path(exists=True)) # need multiple inputs for full flow chart
@click.option("--output", "-o", help="Output filename")
def run_analysis(input, output):
    """Run the biomarker discovery analysis"""
    pass


if __name__ == '__main__':
    run_analysis()