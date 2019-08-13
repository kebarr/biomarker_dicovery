import click # command liine options
# source ~/python3/bin/activate - for python 3
# for now write as command line with input and output names specified

# in input option, nargs is number of input files, we have 4. type tells it its a filepath.
@click.command() 
@click.option("--input", "-i", help="Input files, there must be 4, first is assumed to be the one you want biomarkers for", nargs=4, type=click.Path(exists=True)) # need multiple inputs for full flow chart
@click.option("--output", "-o", help="Output filename")
def run_analysis(input, output):
    """Run the biomarker discovery analysis"""
    pass


if __name__ == '__main__':
    run_analysis()