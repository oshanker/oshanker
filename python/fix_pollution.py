import sys
from pandas import read_csv
from datetime import datetime


# load data
def parse(x):
	return datetime.strptime(x, '%Y %m %d %H')

def main():
    print("-------")
    print(sys.argv)
    # ~/junk/oshanker/oshanker/python/data/
    file_name =  sys.argv[1]     
    dataset = read_csv(file_name+'pollution_raw.csv',  
        parse_dates = [['year', 'month', 'day', 'hour']], index_col=0, date_parser=parse)
    dataset.drop('No', axis=1, inplace=True)
    # manually specify column names
    dataset.columns = ['pollution', 'dew', 'temp', 'press', 'wnd_dir', 'wnd_spd', 'snow', 'rain']
    dataset.index.name = 'date'
    # mark all NA values with 0
    dataset['pollution'].fillna(0, inplace=True)
    # drop the first 24 hours
    dataset = dataset[24:]
    # summarize first 5 rows
    print('head', dataset.head(5))
    print('tail', dataset.tail(5))
    # save to file
    dataset.to_csv(file_name + 'pollution.csv')
    

if __name__   == '__main__':
     main()