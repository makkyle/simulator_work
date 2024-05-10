import pandas as pd
import numpy as np
import os
import sys
import re
import decimal
import subprocess

## import the heavy, light and antigen data 
def extract_data(data_path):
    data = pd.read_excel(os.path.join(data_path,'Nature_LTE.xlsx'),engine='openpyxl')
    data_df = pd.DataFrame(data)

    return data_df

## check if there is binding score for that row of data
def check_data(data,index):
    if  pd.isna(data['binding_score'][index]):
        # type(data['binding_score'][index]) != float and
        return True
    else:
        return False
        
def extract_columns(data,index):
    heavy = data['Heavy'][index]
    light = data['Light'][index]
    antigen = data['variant_seq'][index]
    antigen = antigen.split('\n')
    antigen = antigen[0]
    return [heavy,light,antigen]

def change_path_step1(data_path):
    os.chdir(os.path.join(data_path,'XBCR_net_main'))
    # os.chdir('XBCR-net-main')

heavy = ''
light = ''
antigen = ''

def pass_seq_values(arg1,arg2,arg3):
    global heavy
    global light
    global antigen
    
    heavy = arg1
    light = arg2 
    antigen = arg3

    return heavy,light,antigen

def extract_binding_score(string):
    selected_string = string[-15:]
    potential_float = re.findall("\d+\.\d+",selected_string)
    return potential_float


def process_exponent(value):
    x = str(value)
    ## extract the last two digits of the value (e.g. 10) after the exponent sign (e)
    x = x[-2:]
    x = int(x)
    # print(x)
    # print(f'%.{x}f' % value)
    ## convert the exponent to float number based on the x after e
    return f'%.{x}f' % value

def pass_binding_score(data,index,binding_score):
    data['binding_score'][index] = binding_score
    data.to_excel("C:/Users/harol/Downloads/peter's proposal/YCFS_proposal_material/BCR_cluster/simulator_proof/mutation_pipeline/Nature_LTE.xlsx",engine='openpyxl')
    return data.loc[index,'binding_score'],data

def check_step2(data,index):
    if  pd.isna(data['evolution'][index]):
        # type(data['binding_score'][index]) != float and
        return True
    else:
        return False

def change_path_step2(data_path):
    os.chdir(os.path.join(data_path,'efficient-evolution','efficient-evolution'))

def process_df(df):
    df['freq'] = 0
    for i in range(len(df['output'])):
        df['output'][i] = df['output'][i].replace('#','',1)
        freq = df['output'][i].split(" ")
        print(freq)
        freq_freq = freq[-1]
        print(freq_freq)
        df['freq'][i] = freq_freq
        # if i == 0:
        #     df['mutation'][i] = freq[0]
        # else:
        df['output'][i] = freq[0]
    return df

def convert_str(df):
    obj = df['output'][0]
    return obj

def separateNumbersAlphabets(str):
    numbers = re.findall(r'[0-9]+', str)
    alphabets = re.findall(r'[a-zA-Z]+', str)
    # print(*numbers)
    # print(*alphabets)
    return numbers,alphabets

def build_df(separator_res,df):
    df['original_aa'] = list(np.zeros(int(len(separator_res[0])/2)))
    df['sub_into'] = list(np.zeros(int(len(separator_res[0])/2)))
    df['mutation_pos'] =  list(np.zeros(int(len(separator_res[0])/2)))
    df['freq'] = list(np.zeros(int(len(separator_res[0])/2)))
    
    for i in range(len(separator_res[0])):
        if i%2 ==0:
           df['mutation_pos'][i] = separator_res[0][i]
           df['original_aa'][i] = separator_res[1][i]
        elif i%2 !=0:
           df['freq'][i-1] =  separator_res[0][i]
           df['sub_into'][i-1] = separator_res[1][i]
    ## remove the rows with 0.0 in it 
    df= df[df['original_aa'] != 0]
    df.reset_index(inplace=True)
    df = df.drop(columns=['index'])

    ## check the maximum iteration in "freq" column
    iteration = max(df['freq'])
    return df,iteration

def mutation_generator(org_seq,iteration,df):
    seq = str(org_seq)
    print('in mutation_generator function')
    print(seq)
    print(df)
    origin_seq = seq
    in_loop_count = 0
    mut_count = 0
    mut_seq_list = []
    iteration = int(iteration) + 1
    current_seq = []
    # seq_str = seq
    seq = list(seq)
    # print(seq)

    for i in range(iteration):
        for j in range(len(df)):
        
            if mut_count >0 and mut_count <= int(df['freq'][j]):
                seq_pos = df['mutation_pos'][j]
                seq_pos = int(seq_pos)
                print(seq_pos)
            
                print('seq_before ', seq)
                print(type(seq))
            
                print('before_aa_is :', seq[int(seq_pos)-1])
                seq[seq_pos-1] = df['sub_into'][j]
            
                print('sub into :', seq[seq_pos-1])
                sub_into = seq[seq_pos-1]
            
                print('after_aa_is :', sub_into)
                current_seq = seq
            
                print('seq_after ', current_seq)
                print(type(current_seq))

                diff_list = list(set(current_seq) - set(seq))            
                print('difference ', diff_list)
            
                print('not exceed', j)
                in_loop_count += 1
                print(current_seq)
                # time.sleep(5)
            else:
                print('exceed mut_count', j)
                continue
            if in_loop_count == 0:
                current_seq = []
                print('in_loop_count == 0')
        
        # print(mut_count)
        
        # time.sleep(2)
        mut_count += 1
        print('current_mut_count')
        print(mut_count)
        in_loop_count = 0
        current_seq = ''.join(current_seq)
        mut_seq_list.append(current_seq)
        seq = origin_seq
        seq = list(seq)

    return mut_seq_list


def variantHL_generator_df(df,mut_seq_list,heavy_len):
    ## select particular indices of the list to be saved as csv (for record)
    ## find the largest value in "freq" column of dataframe df 
    max_freq = int(max(df['freq'])) + 1 
    max_freq = int(max_freq)
    range_list = list(range(0,max_freq))
    # print(range_list)

    ## identify the indices that are not in the 'freq' column 
    # freq_df = df['freq'].value_counts
    count_value_list = []
    for i in range(len(df['freq'])):
        # print(df['freq'][i])
        if int(df['freq'][i]) in range_list:
            int_freq = int(df['freq'][i])
            count_value_list.append(int_freq)
        else:
            print("false")

    # print(count_value_list)

    ## remove duplicates from list
    non_dup_list =  list(dict.fromkeys(count_value_list))
    # print(non_dup_list)

    non_dup_list = sorted(non_dup_list)
    print(non_dup_list)
    print(type(non_dup_list))

    len_heavy = heavy_len

    ## save the indices that are useful and not duplicates 
    new_df = pd.DataFrame(columns=['heavy', 'light'])
    new_df['heavy'] = heavy
    new_df['light'] = light

    print(new_df)

    heavy_chain_col = []
    light_chain_col = []

    for i in range(len(non_dup_list)):
        heavy_mut = mut_seq_list[non_dup_list[i]]
        heavy_mut = heavy_mut[:len_heavy]
        print(heavy_mut)
        light_mut =  mut_seq_list[non_dup_list[i]]
        light_mut = light_mut[len_heavy:]
        print(light_mut)
        heavy_chain_col.append(heavy_mut)
        light_chain_col.append(light_mut)

    new_df['heavy'] = heavy_chain_col
    new_df['light'] = light_chain_col
    print(new_df)
    
    return new_df

def binding_score_calculate_s2(H,L,A):
    arg1 = str(H)
    arg2 = str(L)
    arg3 = str(A)
    print('sequences within the function')
    print(H,L,A)
    arg4 = 'XBCR_net'
    arg5 = 'binding'
    arg6 = '0'
    ## run the pred_bcr.py script using subprocess command
    result = subprocess.run(['python','pred_bcr.py','--heavy',arg1,'--light',arg2,'--antig',arg3,
                  '--model_name',arg4,'--data_name',arg5,'--model_num',arg6],capture_output=True,
                          text=True).stdout.strip("\n")
    print('result from subprocess.run')
    print(result)
    print(result[-15:])

    string = result[-15:]
    print('binding_score in string',string)
    binding_score = extract_binding_score(string)[0]
    binding_score = str(binding_score)
    ## convert the exponentials to float numbers
    if "e" in binding_score:
        binding_score = process_exponent(binding_score)
        binding_score = float(binding_score)
    else:
        binding_score = float(binding_score)
    return binding_score
