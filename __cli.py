'''
This is the CLI to access 
'''


import click ##the utility that runs the cli


@click.command()
def start():
    click.echo('Starting message!')

if __name__ == '__main__':
    start()
