import csv 

def read_warmtone(file_location): 
    with open(str(file_location),'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', lineterminator='\r\n')