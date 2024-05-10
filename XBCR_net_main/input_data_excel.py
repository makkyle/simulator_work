import pandas as pd
import numpy as np
import os
import time

from input_data_excel import *

def create_df():
    df = pd.DataFrame()
    return df

def check_csv():
    current_dir = os.getcwd()
    if os.path.exists('binding_score.csv'):
        return True
    else:
        return False

def check_df():
    bs_df = pd.read_csv('binding_score.csv')
    if bs_df['binding_score'].isna() == True:
        return True
    else:
        return False

def input_binding_scores(binding_score):
    
    if check_csv() == True:
        print('the csv file exists')
        if check_df() == True:
            create_df()
            df['binding_score']= [binding_score]
            df.to_csv('binding_score.csv')
        else:
            create_df()
            new_df = pd.DataFrame({'binding_score': [binding_score]})
            df = pd.concat([df,new_df])
            df.to_csv('binding_score.csv')
    elif check_csv() == False:
        print('csv file does not exist')
        df = create_df()
        df['binding_score'] = [binding_score]
        print(df)
        df.to_csv('binding_score.csv')
    time.sleep(10)

    return df
        



