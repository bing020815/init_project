import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import json
import os
from statannotations.Annotator import Annotator
from matplotlib.colors import ListedColormap
from matplotlib_venn import venn3
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import leaves_list
from scipy.spatial.distance import pdist
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.tree import TreeNode
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
from skbio.stats.distance import DistanceMatrix
from itertools import combinations
from matplotlib.lines import Line2D
from scipy.stats import kruskal
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr, spearmanr
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, f_oneway

class PlotBioData:
    """
    Class to generate the statistics of the biological data
    data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
    """
    def __init__(self, data_dict: dict, mapping: dict):
        self.data_dict = data_dict
        self.mapping = mapping
        self.path = './data/plot'
        self.taxa_map = {
            'k': 'Kingdom',
            'p': 'Phylum',
            'c': 'Class',
            'o': 'Order',
            'f': 'Family',
            'g': 'Genus',
            's': 'Species'
        }
        self.sb = StatsBioData(data_dict, mapping)

    def plot_bar(self, df: pd.DataFrame, level: str, target_col: str = '%_hits', title: str = '', show_plot: bool = True):
        """
        function to plot the bar chart of the DataFrame
        df: DataFrame
        level: str, the taxa level
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        title: str, the title of the plot
        """
        # df = df.sort_values(by=target_col, ascending=False)
        # df.plot(kind='bar', x=level, y=target_col, title=title)
        # if show_plot:
        #     plt.show()
        pass

    def plot_stacked(self, groups: list, values_df: pd.DataFrame, title: str = '', show_plot: bool = True):
        """
        #### Values should be expressed as percentage(%)
        function to plot the stacked plot of the DataFrame
        groups: list, the list of the groups
        values_df: DataFrame, the DataFrame with the values
        title: str, the title of the plot
        example:
        groups = ['Group A', 'Group B', 'Group C']
        values_df = pd.DataFrame({
            'otu1': [10, 20, 30],
            'otu2': [20, 30, 40],
            'otu3': [30, 40, 50]
        }, index=['Group A', 'Group B', 'Group C'])
        """
        # preset the plot
        fig, ax = plt.subplots(figsize=(max(10, len(groups) * 0.5), 
                                        max(6, len(values_df.columns) * 0.3)))
        plt.subplots_adjust(right=1) 
        taxa_list = values_df.columns
        filtered_taxa_list = [taxon for taxon in taxa_list if "_" in taxon]
        level = filtered_taxa_list[0].split("_")[0]
        # set the bottom values frame
        bottom_values = [0] * len(groups)

        # check data
        print(f' # taxon: {level}')
        print(f' # groups: {len(groups)}')
        print(f' # values_df: {values_df.shape}')
        output_dir = f'{self.path}/stacked/{level}'
        os.makedirs(output_dir, exist_ok=True)

        # Define a color palette that is colorblind-friendly
        color_palette = [
            '#377eb8',  # Blue
            '#ff7f00',  # Orange
            '#4daf4a',  # Green
            '#f781bf',  # Pink
            '#a65628',  # Brown
            '#984ea3',  # Purple
            '#999999',  # Grey
            '#e41a1c',  # Red
            '#dede00'   # Yellow
        ]

        # Extend the color palette if there are more columns than colors
        if len(values_df.columns) > len(color_palette):
            additional_colors = plt.cm.tab20.colors
            color_palette.extend(additional_colors[:len(values_df.columns) - len(color_palette)])

        for idx, column in enumerate(values_df.columns):
            ax.bar(groups, values_df[column], bottom=bottom_values, label=column, color=color_palette[idx % len(color_palette)])
            bottom_values = [i + j for i, j in zip(bottom_values, values_df[column])]

        # 設定圖例與標題
        ax.set_ylabel(f"Percent of commumity abundance \nat the {self.taxa_map[level]} level(%)")
        ax.set_xlabel('Category')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # tilt the x-axis labels for better readability
        plt.xticks(rotation=45)
        ax.set_title(f'Stacked Plot of {title}')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # 顯示圖表
        # plt.tight_layout()
        # save the plot
        plt.savefig(f'{output_dir}/{title}_{level}_stacked.png', 
                    bbox_inches='tight')
        print(f" Stack plot is saved at {output_dir}/{title}_{level}_stacked.png")
        if show_plot:
            plt.show()
            plt.clf()
            plt.close('all')

    def plot_interactive_stacked(self, groups: list, values_df: pd.DataFrame, title: str = '', show_fig: bool = True):
        """
        #### Values should be expressed as percentage(%)
        function to plot the interactive stacked plot of the DataFrame
        groups: list, the list of the groups
        values_df: DataFrame, the DataFrame with the values
        title: str, the title of the plot

        example:
        groups = ['Category A', 'Category B', 'Category C']

        values_df = pd.DataFrame({
            'otu1': [10, 20, 30],
            'otu2': [20, 30, 40],
            'otu3': [30, 40, 50]
        }, index=['Category A', 'Category B', 'Category C'])
        """
        # preset the plot
        taxa_list = values_df.columns
        filtered_taxa_list = [taxon for taxon in taxa_list if "_" in taxon]
        level = filtered_taxa_list[0].split("_")[0]

        # check data
        print(f' # taxon: {level}')
        print(f' # groups: {len(groups)}')
        print(f' # values_df: {values_df.shape}')
        output_dir = f'{self.path}/stacked/{level}'
        os.makedirs(output_dir, exist_ok=True)

        # transform the DataFrame to long format
        values_df = values_df.reset_index().rename(columns={'index': 'Category'})
        long_df = pd.melt(values_df, id_vars='Category', 
                          var_name='Taxonomy', value_name='Percentage (%)')
        
        # plot the interactive stacked bar plot
        fig = px.bar(long_df, x='Category', y='Percentage (%)', color='Taxonomy', barmode='stack',
                    title=f'Interactive Stacked Bar Plot for {title}',)
        fig.update_layout(barmode='stack', xaxis={'categoryorder':'category ascending'})
        
        # save the plot
        fig.write_html(f"{output_dir}/{title}_{level}_stacked.html")
        print(f" Interactive stacked plot is saved at {output_dir}/{title}_{level}_stacked.html")
        if show_fig: 
            fig.show()
    
    def plot_box(self, values_df: pd.DataFrame, samples_info_df: pd.DataFrame, 
                    plot_single: bool = False, target_x: str = 'Group', 
                    data_type:str = 'otu', show_plot: bool = True,
                    debug_info: bool = False):
        """
        function to plot the boxplot of the DataFrame
        #### boxplot of the data only can be used for the %_hits data
        values_df: DataFrame, the DataFrame with the values
        samples_info_df: DataFrame, the DataFrame with the samples information
        target_x: str, the target column to get the data, default is 'Group'
        example:
        values_df = pd.DataFrame({
            'otu1': [10, 20, 30],
            'otu2': [20, 30, 40],
            'otu3': [30, 40, 50]
        }, index=['sample A', 'sample B', 'sample C'])

        samples_info_df = pd.DataFrame({
            'Group': ['Group A', 'Group B', 'Group A'],
            'Sample_ID': ['sample A', 'sample B', 'sample C'],
            'Category': ['Category A', 'Category B', 'Category A']
        }, index=['sample A', 'sample B', 'sample C'])
        """
        ## safe check
        data_type_list = ['otu', 'diversity', 'pathway']
        assert data_type in data_type_list, f"Invalid data type, please choose from {data_type_list}"
        if data_type == 'otu':
            check_idx = values_df.index.str.endswith(f'{values_df.index.values[0]}')
            values_df_PA = values_df[check_idx]
            check_PA = values_df_PA.sum(axis=1).values[0]
            assert round(check_PA) == 100, 'Value should use %_hits data expressed in percentage(%)'
            
        ## preset the data
        taxa_list = values_df.columns
        filtered_taxa_list = [taxon for taxon in taxa_list if "_" in taxon]
        level = filtered_taxa_list[0].split("_")[0]
        output_dir = f'{self.path}/boxplot/{data_type}/{level}/{target_x}_{level}/'
        os.makedirs(output_dir, exist_ok=True)

        # transform the DataFrame to long format
        values_df = values_df.copy()
        values_df[f'{target_x}'] = samples_info_df[f'{target_x}'].values
        long_df = pd.melt(values_df, id_vars=f'{target_x}', 
                            var_name='Taxonomy', value_name='Value')
        values_df = values_df.drop(columns=[f'{target_x}'])
        num_x_category = long_df[f'{target_x}'].unique()
        print(f" Debug: {len(num_x_category) = }")
        print(f" Debug: {len(sns.color_palette('gray', len(num_x_category))) = }")
        print(f" Debug: {num_x_category = }")
        
        # Define the pairwise comparison
        pairs = list(combinations(set(samples_info_df[f'{target_x}'].values), 2))
        # print(f" Debug: {pairs = }")

        # plot global setting
        sns.set_theme(style="white", font_scale=1.2, rc={"figure.dpi": 300})  # High-quality settings
        
        ## plot the boxplot
        if plot_single:
            # preset the plot
            # plt.figure(figsize=(18, 24))
            num_plots = len(values_df.columns)
            num_cols = 3
            num_rows = (num_plots + num_cols - 1) // num_cols
            plt.figure(figsize=(6 * num_cols, 4 * num_rows))
            # sns.set_theme(style="white", font_scale=1.2, rc={"figure.dpi": 300})
             # High-quality settings
          
            # plot the single boxplot
            for i, taxa in enumerate(values_df.columns, 1):
                # print(f" Debug: {taxa = }")
                plt.subplot((len(values_df.columns) + 2) // 3, 3, i)
                taxa_data = long_df[long_df['Taxonomy'] == taxa]
                ax = sns.boxplot(data=taxa_data, x=f'{target_x}', y='Value',
                                hue=f'{target_x}',
                                palette=sns.color_palette("gray", len(num_x_category)))
                plt.title(f'{taxa}', fontsize=14)
                
                # Define a-xis and y-axis labels
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.xlabel('')
                plt.ylabel('Relative Abundance(%)', fontsize=14)
                plt.title(f'{taxa}', fontsize=16, fontweight='bold')
                # tilt the x-axis labels for better readability
                plt.xticks(rotation=45)

                # Stats annotation for the boxplot
                annotator = Annotator(ax, pairs, data=taxa_data, x=f'{target_x}', y='Value')
                annotator.configure(test='t-test', loc='inside', 
                    text_format='simple',
                    verbose=2,   hide_non_significant=True
                    # , pvalue_thresholds=[[1e-4, "p ≤ 1e-4"], [1e-3, "p ≤ 0.001"], [1e-2, "p ≤ 0.01"], [0.05, "p ≤ 0.05"]]  
                    #, pvalue_thresholds=[[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"]]  
                    )
                annotator.apply_and_annotate()
                # 手動覆蓋標註格式
                for text in ax.texts:  # 遍歷圖形上的所有標註
                    # print(f" Debug: {text.get_text() = }")
                    # print(f" Debug: {text.get_text().split('p') = }")
                    if "p" in text.get_text():
                        try:
                            # 提取 p-value 並進行格式化
                            pval_str = text.get_text().split("p")[1].strip()
                            pval = float(pval_str[1:].strip())
                            if pval <= 0.001:
                                text.set_text(f"*** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.001
                            elif pval <= 0.01:
                                text.set_text(f"** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.01
                            elif pval <= 0.05:
                                text.set_text(f"* $\\it{{p}}$ {pval_str}") # 修改為 p ≤ 0.05
                            else:
                                text.set_text(f"$\\it{{p}}$ = {pval:.2f}")  # 其他情況顯示 p = x.xx
                        except ValueError:
                            continue
            
            # save the plot
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{target_x}_{level}_boxplot.png',
                        dpi=300, bbox_inches='tight')
            print(f" Boxplot for each group is saved at {output_dir}/{target_x}_{level}_boxplot.png")

            # if show_plot:
            #     plt.show()
            #     plt.clf()
            #     plt.close('all')

        else: # plot the boxplot for each target_x
            # plot the boxplot individually
            count = 0
            for taxa in values_df.columns:
                count += 1
                plt.figure(figsize=(6, 4))
                taxa_data = long_df[long_df['Taxonomy'] == taxa]
                ax = sns.boxplot(data=taxa_data, x=f'{target_x}', y='Value', 
                                 hue=f'{target_x}', 
                                 palette=sns.color_palette("gray", len(num_x_category)),
                                 width=0.6, legend=False)
                
                # define the x-axis and y-axis labels
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.xlabel('')
                plt.ylabel('Relative Abundance', fontsize=14)
                plt.title(f'{taxa}', fontsize=16, fontweight='bold')
                plt.xticks(rotation=45)
                
                # Stats annotation for the boxplot
                annotator = Annotator(ax, pairs, data=taxa_data, x=f'{target_x}', y='Value')
                annotator.configure(test='Mann-Whitney', text_format='simple', loc='inside', 
                    verbose=2, hide_non_significant=True
                    # , pvalue_thresholds=[[1e-4, "p ≤ 1e-4"], [1e-3, "p ≤ 0.001"], [1e-2, "p ≤ 0.01"], [0.05, "p ≤ 0.05"]]  
                    # , pvalue_thresholds=[[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"]]
                    )
                annotator.apply_and_annotate()
                # 手動覆蓋標註格式
                for text in ax.texts:  # 遍歷圖形上的所有標註
                    # print(f" Debug: {text.get_text() = }")
                    # print(f" Debug: {text.get_text().split('p') = }")
                    if "p" in text.get_text():
                        try:
                            # 提取 p-value 並進行格式化
                            pval_str = text.get_text().split("p")[1].strip()
                            pval = float(pval_str[1:].strip())
                            if pval <= 0.001:
                                text.set_text(f"*** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.001
                            elif pval <= 0.01:
                                text.set_text(f"** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.01
                            elif pval <= 0.05:
                                text.set_text(f"* $\\it{{p}}$ {pval_str}") # 修改為 p ≤ 0.05
                            else:
                                text.set_text(f"$\\it{{p}}$ = {pval:.2f}")  # 其他情況顯示 p = x.xx
                        except ValueError:
                            continue

                # save the plot
                plt.tight_layout()
                plt.savefig(f"{output_dir}/{taxa}_boxplot.png", 
                            dpi=300, bbox_inches='tight')
                print(f" Boxplot for {taxa} is saved at {output_dir}/{taxa}_boxplot.png")
                
            print(f" # number of boxplot: {count}")

        if show_plot:
            plt.show()
            plt.clf()
            plt.close('all')

    def plot_LDA(self, values_df: pd.DataFrame, samples_info_df: pd.DataFrame, 
                 target_col: str = 'Group', show_plot: bool = True, show_pairwise: bool = True):
        """
        #### Values should be use raletive abundance data expressed in percentage(%)
        function to plot the LDA plot of the DataFrame, filter the sparse features, and transform the data to log10
        values_df: DataFrame, the DataFrame with the LDA values
        target_col: str, the target column to get the data, default is 'Group'
        
        example:
        values_df = pd.DataFrame({
            'otu1': [10, 20, 30],
            'otu2': [20, 30, 40],
            'otu3': [30, 40, 50]
        }, index=['sample A', 'sample B', 'sample C'])
        
        samples_info_df = pd.DataFrame({
            'Group': ['Group A', 'Group B', 'Group A'],
            'Sample_ID': ['sample A', 'sample B', 'sample C'],
            'Category': ['Category A', 'Category B', 'Category A']
        }, index=['sample A', 'sample B', 'sample C'])
        """
        
        ## Preset the data
        # 分離特徵和分組標籤
        X = values_df.copy() # 特徵資料
        y = samples_info_df[[f"{target_col}"]] # 分組標籤
        print(f" Debug: {X.shape = }, {y.shape = }")
        # 計算分組數和特徵數
        num_groups = y[f"{target_col}"].unique() # 分組數量 (list)
        num_features = X.shape[1]  # 特徵數量

        # 動態計算 n_axis/n_components
        n_axis = min(len(num_groups) - 1, num_features)
        print(f" Debug: {n_axis = }")

        # safe check
        output_dir = f'{self.path}/LEfSe/{target_col}'
        os.makedirs(output_dir, exist_ok=True)

        # check the group maximum number
        if len(num_groups) < 2:
            raise ValueError("Too few groups to plot, please add more groups")
        if len(num_groups) == 2:
            var_above_two = False
            test_type = 'Mann-Whitney U'
        else:
            var_above_two = True
            test_type = 'Kruskal-Wallis'
        if len(num_groups) > 4:
            raise ValueError("Too many groups to plot, please down to 4 groups")

        # 過濾稀疏特徵
        threshold = 0.3  # 設定閾值
        X = X.loc[:, (X != 0).mean(axis=0) > threshold]

        # 對特徵資料進行 log 轉換
        X = X.replace(0, 1e-6)
        X = np.log10(X)
        print("X shape:", X.shape)  # 應該是 (samples,fetures)
        print("y shape:", y[f"{target_col}"].shape)  # 應該是一維 (samples,)

        ## LDA Analysis
        lda = LinearDiscriminantAnalysis()
        X_transformed = lda.fit_transform(X, y)
        categories = X.columns
        accuracy = lda.score(X, y)
        print(f"LDA: {X_transformed.shape = }")
        print(f"LDA Accuracy: {accuracy:.2f}")

        # 計算累積解釋比例
        cumulative_variance = np.cumsum(lda.explained_variance_ratio_)
        print(f"Explained variance ratio: {lda.explained_variance_ratio_}")
        print(f"Cumulative explained variance: {cumulative_variance}")

        # 計算 LDA Score（根據特徵權重與輸出貢獻）
        # 取得 LDA 判別軸數量
        lda_scores = lda.coef_.T  # 每個分類單元對所有判別軸的貢獻
        print("lda.coef_ shape:", lda.coef_.shape)
        print("lda_scores shape:", lda_scores.shape)
        print("Number of groups:", len(np.unique(y[f"{target_col}"])))
        print("Groups:", np.unique(y[f"{target_col}"])) 
        # print(f" Debug: \n{lda_scores = }")
        # resolve the lda_scores shape issue, lda_axis may have same num as groups
        if len(num_groups) < lda.coef_.shape[0]:
            n_axis = lda.coef_.shape[0]
        axes_labels = [f'LD{i+1}' for i in range(n_axis)]

        # 動態判斷每個分類單元的特徵性分組
        characteristic_groups = []
        for scores in lda_scores:
            if var_above_two:
                print(f"Debug: ## {lda.classes_ = }")
                # 找到絕對值最大的軸
                print(f"Debug: ## {scores = }")
                max_axis = np.argmax(np.abs(scores))  
                print(f"Debug: ## {max_axis = }, {np.argmax(np.abs(scores)) = }")
                # 根據最大貢獻軸對應分組
                group = lda.classes_[max_axis]
                print(f"Debug: ## {group = }")         
                characteristic_groups.append(group)
            else:
                print(f"Debug: ## {var_above_two = }, {scores = }")
                # 使用 LD1 的分數（scores[0]）判斷
                characteristic_groups.append(lda.classes_[0] if scores[0] < 0 else lda.classes_[1])

        # Significance Test
        p_values_overall = []
        pairwise_significance = []

        for taxa in categories:
            group_data = [X[(y[f"{target_col}"] == group).values][taxa] for group in num_groups]
            # print(f" Debug: ## {taxa = }, {group_data = }")
            # print(f" Debug: ## {type(group_data) = }")
            # print(f" Debug: ## {num_groups = }")   

            # Stats Test 檢定
            try:
                if var_above_two: # 2 組以上
                    # 進行 overall 檢定
                    _, p_overall = kruskal(*group_data)
                    # 進行 pairwise 檢定
                    for group1, group2 in set(combinations(num_groups, 2)):
                        print(f"Debug: ## {taxa = }")
                        print(f"Debug: ## {num_groups = }")
                        print(f"Debug: ## {group1 = }, {group2 = }")
                        data1 = X[(y[f"{target_col}"] == group1).values][taxa]
                        data2 = X[(y[f"{target_col}"] == group2).values][taxa]
                        if data1.equals(data2):
                            print(f"All numbers are identical for {taxa} between {group1} and {group2}, skipping pairwise test.")
                            continue
                        else:
                            _, p_pairwise = kruskal(data1, data2)
                            pairwise_significance.append((taxa, group1, group2, p_pairwise))
                else: # 2 組
                    _, p_overall = mannwhitneyu(*group_data)
                p_values_overall.append(p_overall)
            except ValueError as e:
                print(f"Error in Statistic test for {taxa}: {e}")
                p_values_overall.append(None) 

        # 創建 LDA 結果表
        lda_results = pd.DataFrame({
            'Taxa': X.columns,
            **{label: lda_scores[:, i] for i, label in enumerate(axes_labels)},
            'Characteristic_Group': characteristic_groups,
            'p_value_overall': p_values_overall
        })
        # 計算最大貢獻分數
        lda_results['Max_LDA_Score'] = lda_results[axes_labels].apply(lambda row: max(row, key=abs), axis=1)
        lda_results.to_csv(f'{output_dir}/{target_col}_lda_results.csv', index=False)
        print(f" LDA results saved at {output_dir}/{target_col}_lda_results.csv")

        ## LDA Analysis Plot
        # 動態生成顏色字典
        colors = plt.colormaps['tab10']  # 使用 Matplotlib colormap 生成顏色
        group_colors = {group: colors(i / len(lda.classes_)) for i, group in enumerate(lda.classes_)}

        # 繪製條形圖
        plt.figure(figsize=(12, len(lda_results) * 0.5))
        # bars = plt.barh(lda_results['Taxa'], lda_results['LDA_Score'], color=colors)
        bar_colors = lda_results['Characteristic_Group'].map(group_colors)  # 根據分組動態設置顏色
        bars = plt.barh(lda_results['Taxa'], lda_results['Max_LDA_Score'], color=bar_colors)
        
        # 動態生成 legend
        legend_labels = [plt.Line2D([0], [0], color=color, lw=4, label=group) for group, color in group_colors.items()]
        # plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
        # 動態調整圖例位置，避免與標記重疊
        plt.legend(handles=legend_labels
                   , loc='lower right', bbox_to_anchor=(1.225, -0.5)
                   , fontsize=12
                   , title=f"{target_col}")

        # 添加顯著性標記
        for index, row in lda_results.iterrows():
            p_overall = row['p_value_overall']
            significance = ""
            if p_overall <= 0.001:
                significance = "***"
                p_eval = f"$\\it{{p}}$ ≤ 0.001"
            elif p_overall <= 0.01:
                significance = "**"
                p_eval = f"$\\it{{p}}$ ≤ 0.01"
            elif p_overall <= 0.05:
                significance = "*"
                p_eval = f"$\\it{{p}}$ ≤ 0.05"
            else:
                significance = ""
                p_eval = f"$\\it{{p}}$ = {round(p_overall,3)}"
            
            # 根據 LDA 分數動態調整標記位置
            # 固定文字位置在圖的右側
            text_x = plt.gca().get_xlim()[1] + 0.8  # 獲取當前 x 軸的最大值並偏移 5 單位
            plt.text(text_x, index
                    , f"({p_eval}) {significance}"
                    , va='center'
                    , fontsize=10, color='brown', fontweight='normal')

            # 顯示 pairwise 檢定結果, 標記位置在 0 軸右邊
            pairwise_x = plt.gca().get_xlim()[0] + 5  # 在 0 軸右側偏移 2 單位
            offset = 0  # 偏移量，避免標記重疊
            if len(num_groups) > 2:
                offset = 0  # 偏移量
                for (taxa, group1, group2, p_pairwise) in pairwise_significance:
                    if taxa == row['Taxa']:
                        pair_significance = ""
                        if p_pairwise < 0.001:
                            pair_significance = "***"
                        elif p_pairwise < 0.01:
                            pair_significance = "**"
                        elif p_pairwise < 0.05:
                            pair_significance = "*"
                        if pair_significance and show_pairwise:
                            plt.text(pairwise_x, index - offset
                                , f"{group1} vs {group2}: {pair_significance}"
                                , va='bottom', fontsize=8, color='blue')
                            offset += 0.28  # 增加偏移量以避免重疊

        # 添加標籤與標題
        plt.axvline(0, color='gray', linestyle='--')  # 中心線
        plt.xlabel('LDA Score(log10)', fontsize=12)
        plt.yticks(fontsize=10)
        plt.ylabel('Taxa', fontsize=12)
        plt.title(f"LEfSe LDA for {len(num_groups)} groups [{test_type}]\n")
        plt.tight_layout()

        # 保存圖片
        plt.savefig(f"{output_dir}/{target_col}_lda.png")
        print(f" LDA plot saved at {output_dir}/{target_col}_lda.png")
        if show_plot:
            plt.show()
            plt.clf()
            plt.close('all')

        ## Pairwise Significance           
        if len(num_groups)>2:
            # convert pairwise_significance to DataFrame
            pairwise_sig_df = pd.DataFrame(pairwise_significance, 
                            columns=['Taxa', 'Group1', 'Group2', 'p_value'])
            pairwise_sig_df.to_csv(f'{output_dir}/{target_col}_pairwise_sig.csv', index=False)
            print(f" Pairwise significance saved to {output_dir}/{target_col}_pairwise_sig.csv")
        else:
            print("No Pairwise significance to save")
    
    def plot_annotated_clustermap(self, data_matrix: pd.DataFrame,
                            row_annotations: pd.DataFrame=None, 
                            col_annotations: pd.DataFrame=None,
                            axis:str = "col", alpha: float=0.05, show_plot: bool = True):
        """
        function to plot the Annotated Hierarchical Clustermap of the DataFrame
        data_matrix: DataFrame, the DataFrame with the values
        samples_info_df: DataFrame, the DataFrame with the samples information
        row_annotations: DataFrame, the DataFrame with the row annotations
        col_annotations: DataFrame, the DataFrame with the column annotations
        cmap: str, the color map for the heatmap, default is 'coolwarm'
        figsize: tuple, the size of the figure, default is (12, 8)
        z_score: int, the z-score normalization, default is None
        axis: str, the axis for the cluster, default is 'col'
        alpha: float, the alpha value for the significance, default is 0.05
        """
        ## safe check
        output_dir = f'{self.path}/annotated_clustermap/'
        os.makedirs(output_dir, exist_ok=True)
        # 確保行列標註顏色對應 DataFrame 的順序
        row_colors, col_colors = None, None
        # 生成列分類的顏色映射pool
        if col_annotations is not None:
            col_annotation_colors, used_colors = {}, []
            for col_name in col_annotations.columns:
                colors_dict = self.sb.get_annotation_colors(col_annotations.loc[:,col_name].values,
                                                        used_colors)
                used_colors.extend(colors_dict.values())
                col_annotation_colors[col_name]= colors_dict
                # print(f" Debug: {col_annotation_colors = }")
        if row_annotations is not None:
            # print(f" Debug: {row_annotations}")
            row_annotation_colors = {}
            for row_name in row_annotations.columns:
                colors_dict = self.sb.get_annotation_colors(row_annotations.loc[:,row_name].values)
                row_annotation_colors[row_name]= colors_dict
                # print(f" Debug: {row_annotation_colors = }")
        # 分配行列顏色
        if row_annotations is not None and row_annotation_colors is not None:
            row_colors = pd.DataFrame({
                col: row_annotations[col].map(row_annotation_colors[col])
                for col in row_annotations.columns
            })
        if col_annotations is not None and col_annotation_colors is not None:
            col_colors = pd.DataFrame({
                col: col_annotations[col].map(col_annotation_colors[col])
                for col in col_annotations.columns
            })
        
        # 繪製 Clustermap
        # 重置為默認主題
        sns.reset_defaults()
        #@# dynamic figsize should be adjusted based on the data_matrix shape
        figsize = (data_matrix.shape[0]* 0.8 , 
                   data_matrix.shape[1]* 1.2 if data_matrix.shape[1] < 20 else
                   data_matrix.shape[1]* 4.5)

        g = sns.clustermap(
            data_matrix,
            cmap="coolwarm",
            figsize=(figsize[0] + 2, figsize[1]),  # 增加圖的寬度
            z_score= None,
            row_colors=row_colors,
            col_colors=col_colors,
            dendrogram_ratio=(0.1, 0.1),
            cbar_kws={'label': 'Z-scores'}
        )

        # 計算檢定結果
        ## 進行統計檢定
        significance = self.sb.get_clustermap_test(data_matrix, 
                col_annotations=col_annotations, 
                row_annotations=row_annotations,
                axis=axis, alpha=alpha)
        ## 顯示 ANOVA 檢定結果
        if axis == "row" and significance is not None:
            for row_idx, row_name in enumerate(data_matrix.index):
                g.ax_heatmap.text(
                    -0.5, row_idx + 0.5,  # x 和 y 位置
                    significance.get(row_name, ""),  # 顯示檢定結果
                    ha='center', va='center', fontsize=12, color='black'
                )
        elif axis == "col" and significance is not None:
            for col_idx, col_name in enumerate(data_matrix.columns):
                g.ax_heatmap.text(
                    col_idx + 0.5, -0.5,  # x 和 y 位置
                    significance.get(col_name, ""),  # 顯示檢定結果
                    ha='center', va='center', fontsize=12, color='black'
                )
        # 調整標籤樣式
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

        # 移動 Heatmap 的 Colorbar 到右側
        g.cax.set_position([1.25, 0.7, 0.03, 0.2])  # [x, y, 寬度, 高度]

        # 手動添加 row 顏色圖例（放置在右側）
        if row_annotations is not None and row_annotation_colors is not None:
            row_legend_patches = []
            for annotation, color_mapping in row_annotation_colors.items():
                # 創建 Legend 元素
                row_legend_patches.extend(
                    Patch(facecolor=color, edgecolor='black', label=f"{annotation}: {label}")
                    for label, color in color_mapping.items()
                )
            # 添加分類圖例到右側（與 Heatmap 的 Colorbar 並排）
            g.ax_heatmap.legend(
                handles=row_legend_patches,
                loc='upper left',
                bbox_to_anchor=(1.3, 0.7),  # 調整圖例位置
                fontsize=10,
                title="Row Annotations",
                title_fontsize=12
            )
        # 手動添加 Column 顏色圖例（放置在右側）
        if col_annotations is not None and col_annotation_colors is not None:
            col_legend_patches = []
            for annotation, color_mapping in col_annotation_colors.items():
                # 創建 Legend 元素
                col_legend_patches.extend(
                    Patch(facecolor=color, edgecolor='black', label=f"{annotation}: {label}")
                    for label, color in color_mapping.items()
                )
                plt.gcf().legend(
                handles=col_legend_patches,
                loc="upper right",
                bbox_to_anchor=(1.2, 0.7),  # **右側調整**
                fontsize=10,
                title="Column Annotations",
                title_fontsize=12
                )
        # 保存圖片
        plt.savefig(f"{output_dir}/annotated_clustermap.png", bbox_inches='tight')
        print(f" Clustermap saved at {output_dir}/annotated_clustermap.png")
        if show_plot:
            plt.show()
            plt.clf()
            plt.close('all')
            sns.reset_defaults()

        # save the matrix
        row_order = g.dendrogram_row.reordered_ind
        col_order = g.dendrogram_col.reordered_ind
        data_matrix.iloc[row_order, col_order].to_csv(f"{output_dir}/annotated_clustermap.csv")
        print(f" Data matrix saved at {output_dir}/annotated_clustermap.csv")

    def plot_taxa_heatmap(self, values_df: pd.DataFrame, samples_info_df: pd.DataFrame, 
                     target_col: str = 'Group',
                     threshold: float = 0.3, plot_type :str = "cluster", 
                     show_plot: bool = True):
        """
        function to plot the heatmap of the DataFrame
        values_df: DataFrame, the DataFrame with the values
        title: str, the title of the plot
        threshold: float, the threshold to filter the sparse features, default is 0.3
        plot_type: str, the type of the plot, default is 'cluster', or can be 'heatmap'
        """
        ## safe check for path saving
        output_dir = f'{self.path}/heatmap/{target_col}'
        os.makedirs(output_dir, exist_ok=True)

        ## preset the plot
        df = values_df.copy()
        df.index = samples_info_df[f'{target_col}']
        df = df.groupby(df.index).sum()
        print(f" Debug: {df.index = }")

        # 計算 Z-score 正規化
        # 過濾標準差不為 0 的欄位
        df_filtered = df.loc[:, df.std() != 0]
        # filter out not sigificant data
        df_filtered = df_filtered.loc[:, (df_filtered != 0).mean() > threshold]
        # 按column行方向計算 z-score
        df_normalized = df_filtered.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        df_normalized.to_csv(f'{self.path}/heatmap/{target_col}/{target_col}_norm_heatmap.csv')
        print(f" Normalized data saved to {self.path}/heatmap/{target_col}/{target_col}_norm_heatmap.csv")
        # confirm the z-score data
        # 檢查均值和標準差
        # print("Row-wise mean (should be close to 0):")
        # print(df_normalized.mean(axis=1))  # 如果是 row-wise，axis=1
        # print("Row-wise std (should be close to 1):")
        # print(df_normalized.std(axis=1))
        print("Column-wise mean (should be close to 0):")
        print(df_normalized.mean(axis=0))  # 如果是 Column-wise，axis=0
        print("Column-wise std (should be close to 1):")
        print(df_normalized.std(axis=0))

        # 階層聚類
        linkage_matrix = linkage(pdist(df_normalized.T), method='average')
        
        # save linkage matrix
        pd.DataFrame(linkage_matrix).to_csv(f'{self.path}/heatmap/{target_col}/{target_col}_linkage_matrix.csv')
        print(f" Linkage matrix saved to {self.path}/heatmap/{target_col}/{target_col}_linkage_matrix.csv")

        ## plot the heatmap
        #- plot the heatmap with cluster
        if plot_type == "cluster":
            plt.figure()
            sns.set_theme(style="whitegrid", rc={"figure.dpi": 300})
            g = sns.clustermap( 
                df_normalized, 
                row_cluster=True,  # 聚類行
                col_linkage=linkage_matrix, 
                cmap="coolwarm",
                z_score=None,  # 關閉 Z-score 標準化
                standard_scale=None,  # 不進行標準化
                figsize=( 14 if df_normalized.shape[1] < 10 else 14 + len(df_normalized.columns) * 0.5, 
                         8 if df_normalized.shape[0] < 10 else df_normalized.shape[1] * 0.25),  # 動態調整圖片長寬
                cbar_kws={'label': 'Z-scores'}, # 添加 color bar label
                # cbar_pos=(0.02, 0.7, 0.02, 0.2), # 調整 colorbar 的位置與大小 (x, y, 寬度, 高度)
                dendrogram_ratio=(0.1, 0.1),#(0.02, 0.2),  # 調整樹狀圖與熱圖的比例
                tree_kws={'linewidths': 0.5},  # 樹狀圖線條寬度
            )

            # # 自動調整子圖間距，防止跑位
            # g.fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

            # change legend fontsize
            g.cax.set_ylabel('Z-scores', fontsize=12)
            # 取得 Heatmap 的位置，確保 Legend 與 Heatmap 對齊
            heatmap_position = g.ax_heatmap.get_position()
            cbar_position = g.cax.get_position()
            # **動態調整 Colorbar** # [x, y, 寬度, 高度]
            g.cax.set_position([cbar_position.x0 - 0.05, 
                                cbar_position.y0 -0.35, 
                                0.01, 0.5])
            print(f" Debug: {cbar_position.x0 = }, {cbar_position.y0 = }")
            # 調整 x 軸標籤方向
            plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), 
                     rotation=90, ha='center', fontsize=10)
            # 調整 Y 軸標籤的字體大小
            plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), 
                     rotation=0,  fontsize=10)
            # 調整行和列之間的間距
            g.ax_heatmap.tick_params(axis='x', which='both', length=3)  # 移除 X 軸的刻度線
            g.ax_heatmap.tick_params(axis='y', which='both', length=3)  # 移除 Y 軸的刻度線
            
            # 添加 x 軸標題
            g.ax_heatmap.set_xlabel("Significant Variables", fontsize=14)
            g.ax_heatmap.set_ylabel(f"{target_col}", fontsize=14)

            # 保存圖片
            g.savefig(f"{self.path}/heatmap/{target_col}/{target_col}_cluster.png")
            print(f" Clustermap saved to {self.path}/heatmap/{target_col}/{target_col}_cluster.png")
            
            if show_plot:
                plt.show()
                # 清理殘留的圖形
                plt.clf()  # 清除當前圖形
                plt.close('all')  # 關閉所有圖形
                sns.reset_defaults() # 重置為默認主題
                

        #- plot the heatmap without cluster
        elif plot_type == "heatmap":
            plt.figure()
            # reprocess the data
            # 确保没有 NaN 或 inf 值
            # df_normalized = df_normalized.replace([np.inf, -np.inf], np.nan).dropna()

            # 階層聚類
            linkage_matrix = linkage(pdist(df_normalized.T), method='average')

            # 動態獲取行順序
            if plot_type == "cluster":
                row_order = g.dendrogram_row.reordered_ind
            else:
                linkage_row = linkage(df_normalized, method='average')
                row_order = leaves_list(linkage_row)

            # 重新排列数据
            dendro = dendrogram(linkage_matrix, labels=df_normalized.columns, no_plot=True)
            col_order = dendro['leaves']

            # 重新排列資料
            df_reordered = df_normalized.iloc[row_order, col_order]

            # 繪製熱圖
            plt.figure(figsize=( 14 if df_normalized.shape[1] < 10 else 14 + len(df_normalized.columns) * 0.5, 
                         8 if df_normalized.shape[1] < 10 else 4 + df_normalized.shape[0] * 0.25),)
            
            sns.heatmap(df_reordered, cmap="coolwarm",
                        cbar_kws={'label': 'Z-scores'},  # 添加 color bar label 
                        annot=True, fmt=".2f",  # 顯示數值，並設定格式
                        annot_kws={"size": 8}  # 設定標籤字體大小
                        )

            # 調整 x 軸標籤方向
            plt.xticks(fontsize=10, rotation=90, ha='center')
            plt.yticks(fontsize=10, rotation=0)

            # 添加標題
            plt.title("Heatmap of Significant Variables", fontsize=14)

            # 添加 x 軸標題
            plt.xlabel("Significant Variables", fontsize=14)
            plt.ylabel(f"{target_col}", fontsize=14)

            # 保存圖片
            plt.savefig(f"{self.path}/heatmap/{target_col}/{target_col}_heatmap.png",
                        bbox_inches='tight')
            print(f" Heatmap saved to {self.path}/heatmap/{target_col}/{target_col}_heatmap.png")

            if show_plot:
                plt.show()
                # 清理殘留的圖形
                plt.clf()  # 清除當前圖形
                plt.close('all')  # 關閉所有圖形
                sns.reset_defaults()

        else:
            raise ValueError("Invalid plot type, please choose 'cluster' or 'heatmap'")
        
    def plot_pathway_heatmap(self, bacteria_df: pd.DataFrame, pathway_df: pd.DataFrame, 
                             method: str = 'pearson', show_plot: bool = True):
        """
        function to plot the heatmap of pathway and bacteria
        calculate the peason correlation between the values and plot the heatmap
        if the star if the value <-0.5 or >0.5, show on the heatmap
        values_df: DataFrame, the DataFrame with the values
        threshold: float, the threshold to filter the sparse features, default is 0.3
        """

        ## safe check for path saving
        output_dir = f'{self.path}/heatmap/pathway/{method}'
        os.makedirs(output_dir, exist_ok=True)
        print(f" Debug: {bacteria_df.shape = }, {pathway_df.shape = }")
        correlation_list = ['pearson', 'spearman']
        assert method in correlation_list, f"Invalid method, please choose from {correlation_list}"

        # calculate the pearson correlation, p-value
        corr_matrix = pd.DataFrame(index = bacteria_df.columns, columns = pathway_df.columns)
        pval_matrix = pd.DataFrame(index = bacteria_df.columns, columns = pathway_df.columns)
        for b in bacteria_df.columns:
            for p in pathway_df.columns:
                if method == 'pearson':
                    corr, pval = pearsonr(bacteria_df[b], pathway_df[p])
                elif method == 'spearman':
                    corr, pval = spearmanr(bacteria_df[b], pathway_df[p])
                corr_matrix.loc[b, p] = corr
                pval_matrix.loc[b, p] = pval
        corr_matrix = corr_matrix.astype(float)  # 轉換為浮點型數據
        pval_matrix = pval_matrix.astype(float)  # 轉換為浮點型數據
        corr_matrix.to_csv(f'{output_dir}/pathway_corr.csv')
        print(f" Correlation matrix saved to {output_dir}/{method}_pathway_corr.csv")
        pval_matrix.to_csv(f'{output_dir}/pathway_pval.csv')
        print(f" P-value matrix saved to {output_dir}/{method}_pathway_pval.csv")
        
        
        # 添加星號標記
        significant_mask = pval_matrix < 0.05  # 設定顯著性閾值
        annot = corr_matrix.round(2).astype(str)
        annot[significant_mask] = annot[significant_mask] + " *"
        
        # plot the heatmap for the correlation
        plt.figure(len(pathway_df) * 0.5, len(bacteria_df) * 0.5)
        sns.set_theme(style="whitegrid", font_scale=1, rc={"figure.dpi": 300})
        sns.heatmap(corr_matrix, cmap="coolwarm", 
                    annot=annot, fmt="",
                    # 添加 color bar label
                    cbar_kws={'label': f'Correlation Coefficient({method})'},  
                    annot_kws={"size": 8}  # 設定標籤字體大小
                    )
        
        # 調整 x 軸標籤方向
        plt.xticks(fontsize=10, rotation=90, ha='center')
        plt.yticks(fontsize=10)

        # 添加標題
        plt.title("Heatmap of Pathway Correlation", fontsize=14)

        # 保存圖片
        plt.savefig(f"{output_dir}/{method}_pathway_heatmap.png",
                    bbox_inches='tight')
        print(f" Heatmap saved to {output_dir}/{method}_pathway_heatmap.png")

        if show_plot:
            plt.show()
            # 清理殘留的圖形
            plt.clf()
            plt.close('all')
    
    @classmethod
    def get_box_by_name(cls, tableau_df: pd.DataFrame, samples_info_df: pd.DataFrame,
                        target_x: str, name: str, output_path: str, showfliers=True,
                        show_plot: bool = True):
        """
        function to get the box plot by the name
        example:
        tableau_df = pd.DataFrame({
            'Taxa A': [0.1, 0.2, 0.3],
            'Taxa B': [0.2, 0.3, 0.4],
            'Taxa C': [0.3, 0.4, 0.5]
        }, index=['Sample A', 'Sample B', 'Sample C'])

        samples_info_df = pd.DataFrame({
            'Group': ['Group A', 'Group B', 'Group A'],
            'Sample_ID': ['sample A', 'sample B', 'sample C'],
            'Category': ['Category A', 'Category B', 'Category A']
        }, index=['sample A', 'sample B', 'sample C'])
        """
        ## preset the data
        level = name.split("_")[0]
        output_dir = f'{output_path}/boxplot/{target_x}/{level}/'
        os.makedirs(output_dir, exist_ok=True)

        # transform the DataFrame to long format
        tableau_df = tableau_df[[f'{name}']].copy()
        tableau_df[f'{target_x}'] = samples_info_df[f'{target_x}'].values
        long_df = pd.melt(tableau_df, id_vars=f'{target_x}', 
                            var_name='Taxonomy', value_name='Value')
       
        num_x_category = long_df[f'{target_x}'].unique()
        print(f" Debug: {len(num_x_category) = }")
        print(f" Debug: {len(sns.color_palette('gray', len(num_x_category))) = }")
        print(f" Debug: {num_x_category = }")

        plt.figure(figsize=(6, 4))
        taxa_data = long_df[long_df['Taxonomy'] == name]
        ax = sns.boxplot(data=taxa_data, x=f'{target_x}', y='Value', 
                    hue=f'{target_x}', showfliers=showfliers,
                    palette=sns.color_palette("gray", len(num_x_category)),
                    width=0.6, legend=False)
        
        # define the x-axis and y-axis labels
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # plt.xlabel(f'{target_x}', fontsize=14)
        plt.ylabel('Relative Abundance', fontsize=14)
        plt.title(f'{name}', fontsize=16, fontweight='bold')
        plt.xticks(rotation=45)
        
        # Define the pairwise comparison
        pairs = list(combinations(set(samples_info_df[f'{target_x}'].values), 2))
        
        # Stats annotation for the boxplot
        annotator = Annotator(ax, pairs, data=taxa_data, x=f'{target_x}', y='Value')
        annotator.configure(test='Mann-Whitney', text_format='simple', loc='inside', 
            verbose=2, hide_non_significant=True
            # , pvalue_thresholds=[[1e-4, "p ≤ 1e-4"], [1e-3, "p ≤ 0.001"], [1e-2, "p ≤ 0.01"], [0.05, "p ≤ 0.05"]]  
            # , pvalue_thresholds=[[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"]]
            )
        annotator.apply_and_annotate()
        # 手動覆蓋標註格式
        for text in ax.texts:  # 遍歷圖形上的所有標註
            # print(f" Debug: {text.get_text() = }")
            # print(f" Debug: {text.get_text().split('p') = }")
            if "p" in text.get_text():
                try:
                    # 提取 p-value 並進行格式化
                    pval_str = text.get_text().split("p")[1].strip()
                    pval = float(pval_str[1:].strip())
                    if pval <= 0.001:
                        text.set_text(f"*** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.001
                    elif pval <= 0.01:
                        text.set_text(f"** $\\it{{p}}$ {pval_str}")  # 修改為 p ≤ 0.01
                    elif pval <= 0.05:
                        text.set_text(f"* $\\it{{p}}$ {pval_str}") # 修改為 p ≤ 0.05
                    else:
                        text.set_text(f"$\\it{{p}}$ = {pval:.2f}")  # 其他情況顯示 p = x.xx
                except ValueError:
                    continue

        # save the plot
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{name}_boxplot.png", 
                    dpi=300, bbox_inches='tight')
        print(f" Boxplot for {name} is saved at {output_dir}/{name}_boxplot.png")

        if show_plot:
            plt.show()
            plt.clf()
            plt.close('all')

    @classmethod
    def get_heatmap(cls, data: pd.DataFrame, output_dir: str = None,
                    cmap: str = 'coolwarm', 
                    figsize: tuple = (12, 8), z_score: int = None):
        """
        function to get the heatmap of the data
        data: DataFrame, the DataFrame with the data
        output_dir: str, the directory to save the heatmap, default is None
        cmap: str, the color map for the heatmap, default is 'coolwarm'
        figsize: tuple, the size of the figure, default is (12, 8)
        z_score: int, the z-score normalization, default is None
        """
        # check the data
        if output_dir is None: 
            output_dir = f'{self.path}/heatmap'
        os.makedirs(output_dir, exist_ok=True)

        # plot the heatmap
        sns.set_theme(style="whitegrid")
        g = sns.clustermap(data, cmap=cmap, figsize=figsize, z_score=z_score)
        plt.savefig(f'{output_dir}/heatmap.png', bbox_inches='tight')
        print(f" Heatmap saved to {output_dir}/heatmap.png")
        plt.show()
        plt.clf()
        plt.close('all')
    

class StatsBioData:
    """
    Class to generate the statistics test of the biological data
    data_dict: dict, the data dictionary with key as the sample name, and value as the DataFrame
    """
    def __init__(self, data_dict: dict, mapping: dict):
        self.data_dict = data_dict
        self.mapping = mapping
        self.path = './data/stats'

    def get_stats(self, df: pd.DataFrame, target_col: str = '%_hits'):
        """
        function to get the statistics of the DataFrame
        df: DataFrame
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        """
        stats = df[target_col].describe()
        return stats

    def get_stats_test(self, df1: pd.DataFrame, df2: pd.DataFrame, target_col: str = '%_hits',
                       method: str = 'ttest'):
        """
        function to get the statistics test of the DataFrame
        df1: DataFrame
        df2: DataFrame
        target_col: str, the target column to get the data, default is '%_hits', or can be "num_hits"
        """
        test_list = ['ttest', 'mannwhitneyu']
        assert method in test_list, f"Invalid method, please choose from {test_list}"
        if method == 'ttest':
            # perform the t-test
            t_stat, p_val = ttest_ind(df1[target_col], df2[target_col])
        elif method == 'mannwhitneyu':
            # perform the mannwhitneyu test
            u_stat, p_val = mannwhitneyu(df1[target_col], df2[target_col])    
        return t_stat, p_val
    
    def get_alpha_diversity(self, data: pd.DataFrame, output_dir: str = None):
        """
        function to get the alpha diversity of the data
        data: DataFrame, the DataFrame with the data
        output_dir: str, the directory to save the alpha diversity, default is None
        return: DataFrame, the DataFrame with the alpha diversity
        """
        # ensure the data is in reads format
        from random import randint
        data_sample_value = data.sum(axis=1).to_list()[randint(0, len(data)-1)]
        assert data_sample_value <= 100, f"Invalid data format, data value should be in reads format but got {data_sample_value}"

        # check the data
        if output_dir is None: 
            output_dir = f'{self.path}/alpha_diversity'
        os.makedirs(output_dir, exist_ok=True)

        ids = data.index
        data_df = data.copy()
        # calculate the alpha diversity for the target data
        shannon_indices = alpha_diversity('shannon', data_df.values, ids=ids)
        simpson_indices = alpha_diversity('simpson', data_df.values, ids=ids)
        chao1_indices = alpha_diversity('chao1', data_df.values, ids=ids)
        data_df = data_df.astype(int)
        # 嘗試計算 ACE 指數，如果失敗則使用 Chao1 指數
        try:
            ace_indices = alpha_diversity('ace', data_df.values, ids=ids)
        except ValueError as e:
            print(f"ACE calculation failed: {e}")
            ace_indices = chao1_indices  # 使用 Chao1 指數來代替

            
        # create the dataframe for the alpha diversity
        alpha_df = pd.DataFrame({
            # 'Sample_ID': ids,
            'Shannon': shannon_indices,
            'Simpson': simpson_indices,
            'Chao1': chao1_indices,
            'ACE': ace_indices
        })
        alpha_df.to_csv(f'{output_dir}/alpha_diversity.csv', index=True)
        print(f" Alpha diversity saved to {output_dir}/alpha_diversity.csv")

        return alpha_df
    
    def get_clustermap_test(self, data_matrix: pd.DataFrame, 
                            output_dir: str = None,
                            row_annotations: pd.DataFrame = None,
                            col_annotations: pd.DataFrame = None, 
                            axis: str = "row", 
                            alpha: float = 0.05):
        """
        function to get the clustermap test of the data
         can only be used when row_annotations, row_annotations are setted
        data_matrix: DataFrame, the DataFrame with the data
        output_dir: str, the directory to save the clustermap test, default is None
        row_annotations: DataFrame, the DataFrame with the row annotations, default is None
        col_annotations: DataFrame, the DataFrame with the column annotations, default is None
        axis: str, the axis to perform the test, default is "row", or can be "col"
        alpha: float, the alpha value to perform the test, default is 0.05
        return: DataFrame, the DataFrame with the clustermap test

        example:
        data_matrix = pd.DataFrame({
            'A': [1, 2, 3],
            'B': [4, 5, 6],
            'C': [7, 8, 9]
        }, index=['sample A', 'sample B', 'sample C'])
        row_annotations = pd.DataFrame({
            'Group': ['Category A', 'Category B', 'Category A']
        }, index=['sample A', 'sample B', 'sample C'])
        col_annotations = pd.DataFrame({
            'Level1': ['Firmicutes', 'Proteobacteria', 'Bacteroidetes'],
            'Level2': ['Lactobacillus', 'Pseudomonas', 'Bacteroides']
        }, index=['Bacteria_A', 'Bacteria_B', 'Bacteria_C'])
        """
        from scipy.stats import f_oneway, ttest_ind
        import json

        # check the data
        if output_dir is None: 
            output_dir = f'{self.path}/clustermap_test'
        os.makedirs(output_dir, exist_ok=True)

        significance = {}
        calculation_dict = {}
        if axis == "row":
            # use the top level of the row annotation
            categories = col_annotations.iloc[:, 0].unique() 
            for row_name in data_matrix.index:
                calculation_dict[row_name] = {}
                groups = [data_matrix.loc[row_name, col_annotations[col_annotations.iloc[:, 0] == group].index].dropna().values for group in categories]
                # print(f" Debug: {row_name = }, {categories = }")
                for i, n in enumerate(groups):
                    calculation_dict[row_name][categories[i]] = [arr.tolist() for arr in groups[i]]
                    # print(f" Debug: {categories[i]}: {groups[i]}")
                    continue
                # 過濾掉數量小於 2 的組別，確保每個 Phylum 至少有 2 個數值
                valid_groups = [g for g in groups if len(g) > 1]
                # print(f"{valid_groups =}")
                if len(valid_groups) == 2:
                    f_stat, p_val = ttest_ind(*valid_groups, nan_policy='omit')
                    calculation_dict[row_name]['t-test'] = (float(f_stat), float(p_val))
                    # print(f" Debug: {f_stat = }, {p_val = }")
                    significance[row_name] = "*" if p_val < 0.05 else "n.s."
                elif len(valid_groups) > 2:
                    f_stat, p_val = f_oneway(*valid_groups)
                    calculation_dict[row_name]['oneway_anova'] = {
                        'f_stat': float(f_stat),
                        'p_val': float(p_val)
                    }
                    # print(f" Debug: {f_stat = }, {p_val = }")
                    significance[row_name] = "*" if p_val < 0.05 else "n.s."
                else:
                    # insufficient data, not significant
                    f_stat, p_val = 'i.d.n.s.', 'i.d.n.s.'
                    calculation_dict[row_name]['insufficient_data'] = (f_stat, p_val)
                    significance[col_name] = "i.d.n.s." 
                    print(f"{row_name}: insufficient data for statistical test")
        elif axis == "col":
            # use the top level of the row annotation
            categories = row_annotations.iloc[:, 0].unique()
            # 如果只有 2 個族群，則使用 t-test，否則使用 ANOVA
            test_method = ttest_ind if len(categories) == 2 else f_oneway
            for col_name in data_matrix.columns:
                calculation_dict[col_name] = {}
                groups = [data_matrix.loc[row_annotations[row_annotations.iloc[:, 0] == group].index, col_name].values for group in categories]
                # print(f" Debug: {col_name = }, {categories = }")
                for i, n in enumerate(groups):
                    # print(f" Debug: {categories[i]}: {groups[i]}")
                    calculation_dict[col_name][categories[i]] = [arr.tolist() for arr in groups[i]]
                # 確保所有 group 至少有 2 個樣本，避免錯誤
                if all(len(g) > 1 for g in groups):
                    f_stat, p_val = test_method(*groups)
                    calculation_dict[col_name]['test_method'] = (float(f_stat), float(p_val))
                    # print(f" Debug: {f_stat = }, {p_val = }")
                    significance[col_name] = "*" if p_val < 0.05 else "n.s."
                else:
                    significance[col_name] = "i.d.n.s." # insufficient data, not significant
                    f_stat, p_val = 'i.d.n.s.', 'i.d.n.s.'
                    calculation_dict[col_name]['insufficient_data'] = (f_stat, p_val)
                    print(f"{col_name}: insufficient data for statistical test")
        else:
            print(f"axis must be 'row' or 'col' but got {axis = }.")
            significance = None
        
        # save the calculation dict 
        # dump with the readable json format
        with open(f'{output_dir}/calculation_dict.json', 'w') as f:
            json.dump(calculation_dict, f, indent=4)
        
        return significance
    
    def get_annotation_colors(self, categories: list, used_colors: list = None):
        """
        function to get the annotation colors of the data
        data: DataFrame, the DataFrame with the data
        target_col: str, the target column to get the data, default is 'Group'
        """
        unique_categories = list(set(categories))  # 獲取唯一分類
        num_categories = len(unique_categories)

        # 初始化已使用的顏色列表
        if used_colors is None:
            used_colors = []

        # 取得 "tab20" colormap for 20 unique colors
        # cmap = plt.get_cmap("tab20") 

        # rainbow for 30 unique colors
        cmap1 = plt.get_cmap("tab20b")
        cmap2 = plt.get_cmap("tab20c")
        camp_fourty = ListedColormap(cmap1.colors+cmap2.colors)

        # # rainbow for 50 unique colors
        # cmap = plt.get_cmap("gist_rainbow")

        # 計算已使用顏色的數量
        used_color_count = len(used_colors)

        # 生成顏色，確保新顏色從未使用的部分開始選擇
        new_colors = []
        for i in range(used_color_count, len(camp_fourty.colors)):
            candidate_color = camp_fourty(i)[:3]  # 取得 RGB
            if candidate_color not in used_colors:  # 確保不重複
                new_colors.append(candidate_color + (1.0,))  # 加上 alpha 通道

            if len(new_colors) == num_categories:
                break  # 取得足夠顏色後停止

        # **測試是否有重複**
        duplicated_colors = set(used_colors) & set(new_colors)
        if duplicated_colors:
            print(f"⚠️ 發現重複顏色: {duplicated_colors}")

        # 分配新顏色
        color_mapping = {category: new_colors[i] for i, category in enumerate(unique_categories)}
        
        return color_mapping
    
    def get_sparsed_filtered(self, data: pd.DataFrame, threshold: float = 0.3):
        """
        function to get the sparsed filtered data
        data: DataFrame, the DataFrame with the data
        threshold: float, the threshold to filter the sparse features, default is 0.3
        """
        # filter out the sparse features
        data_filtered = data.loc[:, (data != 0).mean() > threshold]
        filtered_features = list(set(data.columns) - set(data_filtered.columns))
        return data_filtered, filtered_features
    
    def get_beta_diversity_test(self, data: pd.DataFrame, method: str = 'braycurtis'):
        """
        function to get the beta diversity test of the data
        data: DataFrame, the DataFrame with the data
        method: str, the method to calculate the beta diversity, default is 'braycurtis'
        """
        beta_diversity = beta_diversity(method, data.values)
        return beta_diversity
    def get_filtered_features(self, data_df:pd.DataFrame, samples_info_df:pd.DataFrame,
                          control_group:str,
                          treatment_group:str, 
                          target_col:str, 
                          filter_method:str = 'stats', 
                          note:str = 'pathway'):
        """
        get the filtered features based on the filter method
        return the significant features in list

        filter_method: 'stats', 'variance', 'stats_variance', 'biological'
        note: custom note, can be 'pathway', 'metabolite', 'taxa'...
        """
        assert target_col in samples_info_df.columns, f"{target_col} not in {samples_info_df.columns}"
        assert control_group in samples_info_df[target_col].values, f"{control_group} not in {samples_info_df[target_col].values}"
        assert treatment_group in samples_info_df[target_col].values, f"{treatment_group} not in {samples_info_df[target_col].values}"
        ## safe check for path saving
        output_dir = f'{self.path}/filtered_features'
        os.makedirs(output_dir, exist_ok=True)

        # 統計檢定 (ANOVA 或 t-test),  變異度篩選 (Variance Thresholding),
        # 統計 + 變異度雙重篩選, 生物學意義篩選 (Pathway Enrichment Analysis)
        if filter_method == 'stats':
            p_values = {}
            samples_info_df_copy = samples_info_df.copy().set_index('Sample_ID')
            for pathway in data_df.columns:
                group_1 = data_df.loc[samples_info_df_copy[f'{target_col}'] == f'{control_group}'][pathway]
                group_2 = data_df.loc[samples_info_df_copy[f'{target_col}'] == f'{treatment_group}'][pathway]
                
                # 如果只有兩組，使用 t-test
                stat, p_val = ttest_ind(group_1, group_2)
                p_values[pathway] = {'stats':stat,
                                    'p-value':p_val}
            # 篩選顯著 pathway (p < 0.05)
            significant_featrues = [k for k, v in p_values.items() if v['p-value'] < 0.05]

        num_of_features = len(significant_featrues)
        total_features = len(data_df.columns)
        p_values['stats']={
            'total features': total_features,
            'significant features': num_of_features
            }
        # output json file
        with open(f'{output_dir}/{filter_method}_features_{num_of_features}_{note}.json', 'w') as f:
            json.dump(p_values, f, indent=4)
        print(f"Filtered features is saved in {output_dir}/{filter_method}_features_{num_of_features}_{note}.json")
        print(f"## significant features: {significant_featrues = }")
        print(f"## significant features: {len(significant_featrues) = }")

        return significant_featrues