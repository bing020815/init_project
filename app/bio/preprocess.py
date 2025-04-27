import pandas as pd
from functools import reduce
import os
import glob

class PrepBioData:
    """
    Class to prepare the biological data for visualization
    data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
    mapping: dict, the mapping dictionary with key as the taxa level, and value as the prefix
    """
    def __init__(self, data_dict: dict, mapping: dict):
        self.mapping = mapping
        self.data_dict = data_dict
        self.proc_level = ''
        self.path = os.path.join('data', 'PrepBioData')

    def get_taxa_level(self, df: pd.DataFrame, level: str, target_col: str = '%_hits'):
        """
        function to get the specific taxa level data, and sum the data with the same taxa level
        if the level is phylum, the value will be set with prefix 'p'
        if the level is empty, the value will be set as 'level_Others'
        df: DataFrame
        level: str, the taxa level
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        example:
         level: Phylum
         target_col: '%_hits'
        output:
         Phylum  %_hits
        """
        assert level in ["Kindom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
                         ,"Species_only"], f"Level {level} not in the DataFrame columns"
        cols = [f'{level}', f'{target_col}']
        self.proc_level = level
        # print(f'get_taxa_level: {level}')
        if level in ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']:
            df = df[cols].copy()
            # set the value with prefix 'p' if level is phylum
            df.loc[:, f'{level}']  = f"{self.mapping['taxa'][level]}_" + df[f'{level}' + "_spp."]
            df = df.fillna(f"{self.mapping['taxa'][level]}_Others")
            df = df.groupby(f'{level}').sum().reset_index()
            # print(f" Debug: {df.columns = }")
            df = df[[f'{level}', f'{target_col}']]
        elif level == 'Species':
            # get the Species as s_Genus_Species
            df = df[[f'Genus', f'{level}', f'{target_col}']].copy()
            df['Species'] = df['Genus'] + '_' + df['Species']
            df['Species'] = df['Species'].fillna('Others')
            df['Species'] = df['Species'].map(lambda x: f"s_{x}")
            df = df[cols]
            df = df.groupby(f'{level}').sum().reset_index()
        elif level == 'Species_without_Genus_combination':
            # get the Species as s_Species
            df = df[[f'Species', f'{target_col}']].copy()
            df['Species'] = df['Species'].fillna('Others')
            cols = [f'Species', f'{target_col}']
            df = df[cols]
            df = df.groupby(f'Species').sum().reset_index()
        # check if the values sum to 1
        if target_col == '%_hits':
            assert round(df['%_hits'].sum()) == 1, f"{target_col} values do not sum to 1: {df['%_hits'].sum()}"
            pass
        # print(f" Debug: {df.columns = }")
        df = df.sort_values(by=f'{target_col}', ascending=False)
        df = df.reset_index(drop=True)
        return df

    def get_tabular_data(self, df: pd.DataFrame, level: str, target_col: str = '%_hits'):
        """
        #### Values perform better in "num_hits" than "%_hits" 
        df structure:
        columns: ['p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others']
        index: ['Group1', 'Group2', 'Group3']
        values: percentage of each taxa
        target_col: the target column to get the data, default is '%_hits', or can be "num_hits"

        example:
        df = pd.DataFrame(
                {'p_bacteria1': [50, 40, 30],
                'p_bacteria2': [20, 30, 40],
                'p_bacteria3': [10, 15, 10],
                'Others': [20, 15, 20]
                    }, 
            index=['Group1', 'Group2', 'Group3']
            )
        output:
        if target_col = '%_hits' returns the percentage of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50.0          20.0          10.0          20.0

        if target_col = 'num_hits' returns the number of hits of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50            20            10            20
        """
        # get the data of the specific level
        taxa_df = self.get_taxa_level(df, level=level, target_col=target_col)
        # Transpose the DataFrame
        transposed_df = taxa_df.transpose()
        # Set the first row as the column headers
        transposed_df.columns = transposed_df.iloc[0]
        transposed_df = transposed_df[1:]
        # reset column index name
        transposed_df.columns.name = None
        # reset the index
        transposed_df.reset_index(drop=True, inplace=True)
        if target_col == '%_hits':
            # multiply 100 to get the percentage
            transposed_df = transposed_df * 100
        return transposed_df
    
    def get_all_tabular_data(self, level: str, target_col: str = '%_hits'):
        """
        data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
        level: str, the taxa level
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        
        output:
        if target_col = '%_hits' returns the percentage of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50.0          20.0          10.0          20.0
        
        if target_col = 'num_hits' returns the number of hits of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50            20            10            20
        """
        # get the data for each sample
        values_dict = {}
        count = 0
        index_list = []
        for key, value in self.data_dict.items():
            count += 1
            df = self.get_tabular_data(value, level=level, target_col=target_col)
            values_dict[f'{key}'] = df
            index_list.append(f'{key}')
        print(f' # of samples: {count} samples')
        # concat the dataframes of each sample
        samples_concat = pd.concat(values_dict.values())
        samples_concat.index = index_list
        # Replace NaN values with 0 using map
        samples_concat = samples_concat.map(lambda x: 0 if pd.isna(x) else x)
        print(f' # of bacterium: {samples_concat.shape[1]} bacterium')
        if not os.path.exists(f'./data/PrepBioData'):
            os.makedirs(f'./data/PrepBioData')
        samples_concat.to_csv(f'{self.path}/{level}_{target_col}_tabular_data.csv')
        print(f" Tabular data is saved to {self.path}/{level}_{target_col}_tabular_data.csv")
        return samples_concat
    
    def save_all_tabular_parquet_data(self, data_dict:dict, 
                                      level: str = 'Species_without_Genus_combination', 
                                      target_col: str = '%_hits'):
        """
        data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
        level: str, the taxa level
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        
        output:
        if target_col = '%_hits' returns the percentage of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50.0          20.0          10.0          20.0
        
        if target_col = 'num_hits' returns the number of hits of each taxa
                'p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'p_Others'
         Values      50            20            10            20
        """
        # 巨量資料處理
        parquet_path = f'{self.path}/parquet_files'
        os.makedirs(parquet_path, exist_ok=True)
        # clear the folder before saving the data
        files = glob.glob(f'{parquet_path}/*')
        for f in files:
            os.remove(f)
        count = 0
        # get the data for each sample and save as parquet
        for key, value in data_dict.items():
            count += 1
            df = self.get_tabular_data(value, level=level, target_col=target_col)
            df.fillna(0, inplace=True)
            df.to_parquet(f'{parquet_path}/{key}.parquet')
            print(f" {count}. Tabular parquet data is saved to {parquet_path}/{key}.parquet")
        print(f" All tabular parquet data is saved to {parquet_path}")
        return parquet_path

    def get_all_tabular_data_from_parquet(self, parquet_path, batch_size=100):
        # 動態讀取所有 Parquet 檔案
        parquet_files = glob.glob(f"{parquet_path}/*.parquet")
        print(f" # of parquet files: {len(parquet_files)}")
        counter = 0
        all_batches = []
        index_list = []
        for i in range(0, len(parquet_files), batch_size):
            batch_files = parquet_files[i:i + batch_size]
            df_list = [pd.read_parquet(file) for file in batch_files]
            samples_concat = pd.concat(df_list, keys=[file.split("/")[-1].replace(".parquet", "") for file in batch_files])
            all_batches.append(samples_concat)
            if len(parquet_files) > 100:
                counter += batch_size
            # 保存索引
            for file in batch_files:
                index_list.append(file.split("/")[-1].replace(".parquet", ""))
            print(f" Processed up to {counter} cumulative samples...")

        print(f" { len(all_batches) = }")
        
        # 合併所有批次的資料
        if len(all_batches)>=1:
            samples_concat = pd.concat(all_batches)
            samples_concat.index = index_list
            # fill NaN with 0
            samples_concat.fillna(0, inplace=True)
            print(f" samples_concat shape: {samples_concat.shape}")
            print(f" samples concatenation done")
        else:
            raise ValueError('No parquet files found to be concatenated')
        return samples_concat
    
    def get_all_metabolites(self, metab_dict: dict, value_type: str = '%_hits'):
        """
        Get the metabolite data for the column annotations
        It gets the metabolite data for each sample and merge them into a single DataFrame
        metab_dict: dict, the metabolite data dictionary with key as the sample name, and value as the DataFrame
        value_type: str, the value type to get the data, default is '%_hits', or can be "num_hits"

        output:
        if value_type = 'num_hits' returns the percentage of each metabolite

        if value_type = '%_hits' returns the percentage of each metabolite

        """
        # preset the data
        output_dir = f'{self.path}'
        os.makedirs(output_dir, exist_ok=True)

        # get the data for each sample
        # outer join the metabolite data on pathway as key
        metab_data_cols = len(set(list(metab_dict['batch0'].columns)+
                                list(metab_dict['batch1'].columns)+
                                list(metab_dict['batch2'].columns)))
        metab_data_rows = len(set(list(metab_dict['batch0'].pathway.values)+
                                list(metab_dict['batch1'].pathway.values)+
                                list(metab_dict['batch2'].pathway.values)))
        print(f" # of columns for Metabolite data: {metab_data_cols-2} samples")
        print(f" # of rows for Metabolite data: {metab_data_rows} pathways")
        # 取出所有 dataframe
        dfs = list(metab_dict.values())
        # 利用 reduce 動態依序以 "pathway" 欄位外部合併
        merged_df = reduce(lambda left, right: pd.merge(left, right, on=['pathway', 'description'], how="outer"), dfs)
        assert merged_df.shape[0] == metab_data_rows, f'metabolite data should have the rows of {metab_data_rows}, but got {merged_df.shape[0]}'
        assert merged_df.shape[1] == metab_data_cols, f'metabolite data should have the columns of {metab_data_cols}, but got {merged_df.shape[1]}'
        
        # index setting
        # print(f" {merged_df.columns = }")
        info_df = merged_df[['pathway', 'description']].copy()
        merged_df.index = info_df['description'].values

        if value_type == '%_hits':
            # 提取數值欄位（即樣本數據）
            metab_df = merged_df.drop(columns=['pathway', 'description']).fillna(0).copy()
            # 以樣本為基準計算相對豐度（轉換為百分比）
            # 計算每個樣本的總和，並將每個路徑的數值除以該樣本的總和，再乘以 100 轉為百分比
            metab_df = (metab_df.div(metab_df.sum(axis=0), axis=1) * 100)
            # verify the data, each column should sum up to 100
            assert (round(metab_df.iloc[:, 2:].sum(axis=0)) == 100).all(), 'each column should sum up to 100'
            # transpose the DataFrame
            metab_df = metab_df.T
        elif value_type == 'num_hits':
            metab_df = merged_df.drop(columns=['pathway', 'description']).fillna(0).copy().T

        metab_df.to_csv(f'{output_dir}/{value_type}_metabolite_data.csv', index=True)
        print(f" Metabolite data is saved to {output_dir}/{value_type}_metabolite_data.csv")
        info_df.to_csv(f'{output_dir}/metabolite_info.csv', index=True)
        print(f" Metabolite info is saved to {output_dir}/metabolite_info.csv")

        return metab_df
    
    def get_annotated_taxa(self, keys_of_batches: list):
        """
        Get the annotated taxa DataFrame for the column annotations
        It gets the Phylum, Genus, and Species columns from the data_dict
        It is used to get the column annotations for the Hierarchical annotated cluster heatmap
        """
        # preset the data
        output_dir = f'{self.path}/'
        os.makedirs(output_dir, exist_ok=True)
        key_dict ={}
        for key in keys_of_batches:
            annotated_taxa_df_tmp = self.data_dict[key].loc[:,['Phylum', 'Genus', 'Species']].copy()
            annotated_taxa_df_tmp['Species'] = annotated_taxa_df_tmp['Genus'] + '_' + annotated_taxa_df_tmp['Species']
            annotated_taxa_df_tmp['Species'] = annotated_taxa_df_tmp['Species'].fillna('Others')
            annotated_taxa_df_tmp['Species'] = annotated_taxa_df_tmp['Species'].map(lambda x: f"s_{x}")
            annotated_taxa_df_tmp.index = annotated_taxa_df_tmp['Species'].values
            key_dict[key] = annotated_taxa_df_tmp
        # concat the dataframes of each sample from dict
        annotated_taxa_df = pd.concat(key_dict.values())
        annotated_taxa_df.to_csv(f'{self.path}/col_annotations.csv')
        annotated_taxa_df = annotated_taxa_df.fillna('Others').drop(columns=['Species'])
        # # get the annotated taxa
        # annotated_taxa_df = self.data_dict[first_key].loc[:,['Phylum', 'Genus', 'Species']].copy()
        # annotated_taxa_df['Species'] = annotated_taxa_df['Genus'] + '_' + annotated_taxa_df['Species']
        # annotated_taxa_df['Species'] = annotated_taxa_df['Species'].fillna('Others')
        # annotated_taxa_df['Species'] = annotated_taxa_df['Species'].map(lambda x: f"s_{x}")
        # annotated_taxa_df.index = annotated_taxa_df['Species'].values
        # annotated_taxa_df.to_csv(f'{self.path}/col_annotations.csv')
        # annotated_taxa_df = annotated_taxa_df.fillna('Others').drop(columns=['Species'])
        # print(f" Annotated taxa saved to {self.path}/col_annotations.csv")
        return annotated_taxa_df

# def get_taxa_level(df: pd.DataFrame, level: str, target_col: str = '%_hits', mapping: dict = mapping):
#     """
#     function to get the specific taxa level data, and sum the data with the same taxa level
#     if the level is phylum, the value will be set with prefix 'p'
#     if the level is empty, the value will be set as 'Others'
#     df: DataFrame
#     level: str, the taxa level
#     target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
#     """
#     cols = [f'{level}', '%_hits']
#     # print(f'get_taxa_level: {level}')
#     df = df[cols]
#     # set the value with prefix 'p' if level is phylum
#     df.loc[:, f'{level}']  = f"{mapping['taxa'][level]}_" + df[f'{level}']
#     df = df.fillna('Others').groupby(f'{level}').sum().reset_index()
#     df = df[[f'{level}', f'{target_col}']]
#     df = df.sort_values(by=f'{target_col}', ascending=False)
#     df = df.reset_index(drop=True)
#     return df


# def get_tabular_data(df: pd.DataFrame, level: str, target_col: str = '%_hits'):
#     """
#     df structure:
#      columns: ['p_bacteria1', 'p_bacteria2', 'p_bacteria3', 'Others']
#      index: ['Group1', 'Group2', 'Group3']
#      values: percentage of each taxa
#      target_col: the target column to get the data, default is '%_hits', or can be "num_hits"

#     example:
#     df = pd.DataFrame(
#             {'p_bacteria1': [50, 40, 30],
#              'p_bacteria2': [20, 30, 40],
#              'p_bacteria3': [10, 15, 10],
#              'Others': [20, 15, 20]
#                 }, 
#         index=['Group1', 'Group2', 'Group3']
#         )
#     """
#     # get the data of the specific level
#     taxa_df = get_taxa_level(df, level=level)
#     # Transpose the DataFrame
#     transposed_df = taxa_df.transpose()
#     # Set the first row as the column headers
#     transposed_df.columns = transposed_df.iloc[0]
#     transposed_df = transposed_df[1:]
#     # reset column index name
#     transposed_df.columns.name = None
#     # reset the index
#     transposed_df.reset_index(drop=True, inplace=True)
#     # multiply 100 to get the percentage
#     transposed_df = transposed_df * 100
#     return transposed_df

# def get_all_tabular_data(data_dict: dict, level: str, target_col: str = '%_hits'):
#     """
#     data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
#     level: str, the taxa level
#     target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
#     """
#     # get the data for each sample
#     values_dict = {}
#     count = 0
#     index_list = []
#     for key, value in data_dict.items():
#         count += 1
#         df = get_tabular_data(value, level=level, target_col=target_col)
#         values_dict[f'{key}'] = df
#         index_list.append(f'{key}')
#     print(f'# of samples: {count} samples')
#     # concat the dataframes of each sample
#     samples_concat = pd.concat(values_dict.values())
#     samples_concat.index = index_list
#     samples_concat = samples_concat.fillna(0)
#     return samples_concat