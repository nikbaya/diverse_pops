#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 10:17:00 2020

@author: nbaya
"""

import hail as hl
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

plots_dir = '/Users/nbaya/Downloads'
data_dir = '/Users/nbaya/Downloads'

clump_pops_dict = {'100001.prs.tsv.gz':['AMR', 'CSA', 'EAS', 'EUR', 'MID'],
                   '30600.prs.tsv.gz':['AFR', 'AMR', 'EAS', 'EUR', 'MID'],
                   '100320.prs.tsv.gz':['EUR'],
                   '100020.prs.tsv.gz':['AFR','AMR', 'EAS', 'EUR', 'MID'],
                   '250.2.prs.tsv.gz':['AFR', 'CSA', 'EAS', 'EUR', 'MID'],
                   'E78.prs.tsv.gz':['AFR', 'AMR', 'EAS', 'EUR', 'MID'],
                   }    

def plot_prs_corr_hist(pop=None, logscale=False):

    df0 = pd.read_csv(f'{data_dir}/pilot_assess_prs{"" if pop is None else f".{pop}" }.tsv.gz', 
                     sep='\t', compression='gzip')
    df1 = df0[~df0.prs_corr_r2.isna()]
    trait_types = sorted(set(df1.trait_type))
    
    bins = np.linspace(0,1,51)
    x = bins[:-1]+(bins[1:]-bins[:-1])/2
    width = (bins[1:]-bins[:-1])[0]
    prev = 0
    for trait_type in trait_types:
        df = df1[df1.trait_type==trait_type]
        n_curr, _ = np.histogram(df.prs_corr_r2)
        plt.bar(x=)
    
    plt.figure(figsize=(6*1.5, 4*1.5))
    plt.hist(df[~df.prs_corr_r2.isna()].prs_corr_r2, bins=np.linspace(0,1,101))
    plt.xlim([0,1])
    plt.xlabel(r'PRS $r^2$')
    plt.ylabel('density')
    if logscale:
        plt.yscale('symlog', linthreshy=4)
    plt.title(f'PRS r2 for {} pop-pheno-pval combinations'+
              f'{"" if pop is None else f" (pop: {pop})"}\n'+
              f'({df[df.prs_corr_r2.isna()].shape[0]} with r2=NA excluded)')
    plt.savefig(f'{plots_dir}/prs_corr_hist{"" if pop is None else f".{pop}" }{".logscale" if logscale else ""}.png', dpi=300)

def plot_prs_corr_holdout(pop):
    df0 = pd.read_csv(f'{data_dir}/pilot_assess_prs{"" if pop is None else f".{pop}" }.tsv.gz', 
                     sep='\t', compression='gzip')
    df1 = df0[(~df0.prs_corr_r2.isna())&(~df0.clump_pops_str.str.contains(pop))]
    df2 = df1.groupby('p_threshold')['prs_corr_r2','prs_corr_pval'].mean()
    df3 = df1.groupby('p_threshold')['prs_corr_r2','prs_corr_pval'].sem()
#    df4 = df1.groupby('p_threshold')['prs_corr_r2','prs_corr_pval'].std()
    df5 = df1.groupby('p_threshold')['prs_corr_r2'].count()

    fig, (ax0, ax1) = plt.subplots(1,2,figsize=(12*1.2,4*1.5))
    for metric, ax in [('r2', ax0), ('pval',ax1)]:
        ax.errorbar(x=df2.index, y=df2[f'prs_corr_{metric}'], 
                    yerr=2*df3[f'prs_corr_{metric}'], fmt='.-')
        if metric=='r2':
            ax.set_ylabel(r'$r^2$')
        elif metric=='pval':
            ax.set_ylabel(r'$r^2$ pval')
            ax.set_ylim([0,1])
        ax.set_xscale('symlog',linthreshx=1e-9)
        ax.set_xlabel('clumping pval threshold')
        ax.set_xlim(ax.get_xlim()[::-1])
    plt.suptitle(f'Holdout population {pop}\np_threshold phen cts: {dict(df5[::-1])} ')
    plt.savefig(f'{plots_dir}/prs_corr_holdout.{pop}.png',dpi=300)
    
def plot_score_phen(fname, force_quant=False, logscale=False, pop=None, p_threshold=None):
    r'''
    For plotting scatter plot of phenotype vs. score (for quantitative traits 
    or if force_quant=True) or a histogram stratified by case/control status
    for binary phenotypes
    '''
    df0 = pd.read_csv(f'{data_dir}/{fname}', sep='\t', compression='gzip')
    trait_type, phenocode, phen_desc = df0[['trait_type','phenocode','description']].values[0]
    phen_desc = ' '.join(phen_desc.split(' ')[:4])+('...' if len(phen_desc.split(' '))>4 else '')
    df1 = df0[~df0.both_sexes.isna()]
    pops = sorted(set(df1['pop'])) if pop is None else [pop]
    clump_pops = clump_pops_dict[fname]
    r2_all = np.corrcoef(df1.both_sexes, df1.score)[0,1]**2
    p_thresholds = sorted(set(df1.p_threshold)) if p_threshold is None else [p_threshold]
    for p_threshold in p_thresholds:
        df2 = df1[df1['p_threshold']==p_threshold]
        for pop in pops:
            df = df2[df2['pop']==pop]
            r_pop, pval_pop = stats.pearsonr(df.both_sexes, df.score)
            r2_pop = r_pop**2
            plt.figure(figsize=(6*1.5, 4*1.5))
            if trait_type in ['continuous','biomarkers'] or force_quant:
                plt.plot(df.both_sexes, df.score, '.', alpha=0.5)
                m, b = np.polyfit(df.both_sexes, df.score, 1)
                plt.plot(df.both_sexes, m*df.both_sexes+b, 'k--')
                plt.xlabel('phenotype')
                plt.ylabel('score')
            else:
                min_score = df.score.min()
                max_score = df.score.max()
                bins = np.linspace(min_score, max_score, 100)
                controls_hist, _ = np.histogram(df[df.both_sexes==0].score, bins=bins)
                cases_hist, _ = np.histogram(df[df.both_sexes==1].score, bins=bins)
                x = bins[:-1]+(bins[1:]-bins[:-1])/2
                width = (bins[1:]-bins[:-1])[0]
        
                plt.bar(x=x, height=controls_hist, width=width)
                plt.bar(x, cases_hist, bottom=controls_hist, width=width)
                if logscale:# or min(max(controls_hist), max(cases_hist))<max(max(controls_hist), max(cases_hist))/100:
                    plt.yscale('symlog', linthreshy=5)
                plt.xlabel('score')
                plt.ylabel('count')
                plt.legend(['controls','cases'])
            plt.title(f'{phen_desc} (phenocode:{phenocode}, clump pops: {",".join(clump_pops)}, pval={p_threshold})\n'+
                      f'{pop} r2={round(r2_pop,3)}, pval={round(pval_pop,3)} (all pop r2={round(r2_all,3)}), {df.shape[0]} samples '+
                      f'({df0[df0.both_sexes.isna()].shape[0]} missing pheno data)')
            plt.savefig(f'{plots_dir}/score_pheno_plot.{trait_type}.{phenocode}.{pop}.pval_{p_threshold}'+
                        f'{".logscale" if logscale else ""}{".force_quant" if force_quant else ""}.png',dpi=300)    

def plot_r2_pval_threshold(fname, pop):
    df0 = pd.read_csv(f'{data_dir}/{fname}', sep='\t', compression='gzip')
    trait_type, phenocode, phen_desc = df0[['trait_type','phenocode','description']].values[0]
    phen_desc = ' '.join(phen_desc.split(' ')[:4])+('...' if len(phen_desc.split(' '))>4 else '')
    df1 = df0[~df0.both_sexes.isna()]
    clump_pops = clump_pops_dict[fname]
    p_thresholds = sorted(set(df1.p_threshold))
    r2_pval_list = []
    for p_threshold in p_thresholds:
        df = df1[(df1['p_threshold']==p_threshold)&(df1['pop']==pop)]
        r, r_pval= stats.pearsonr(df.both_sexes, df.score)
        r2 = r**2
        r2_pval_list.append((r2, r_pval))

    fig, (ax0, ax1) = plt.subplots(1,2,figsize=(12*1.2,4*1.5))
    for metric, ax in [('r2', ax0), ('r2_pval',ax1)]:
        if metric=='r2':
            ax.plot(p_thresholds, [r2 for r2,pval in r2_pval_list], '.-')
            ax.set_ylabel(r'$r^2$')
        elif metric=='r2_pval':
            ax.plot(p_thresholds, [pval for r2,pval in r2_pval_list], '.-')
            ax.set_ylabel(r'$r^2$ pval')
            ax.set_ylim([0,1])
        ax.set_xscale('symlog',linthreshx=1e-9)
        ax.set_xlabel('clumping pval threshold')
        ax.set_xlim(ax.get_xlim()[::-1])
    plt.suptitle(f'{phen_desc} (phenocode:{phenocode}, clump pops: {",".join(clump_pops)})\n'+
                      f'{pop} , {df.shape[0]} samples '+
                      f'({df0[df0.both_sexes.isna()].shape[0]} missing pheno data)')
    plt.savefig(f'{plots_dir}/r2_pval_threshold.{trait_type}.{phenocode}.{pop}.png',dpi=300)    
            
def test_prs_score_diff(fname):
    df0 = pd.read_csv(f'{data_dir}/{fname}', sep='\t', compression='gzip')
    trait_type, phenocode, phen_desc = df0[['trait_type','phenocode','description']].values[0]
    df1 = df0[~df0.both_sexes.isna()]
    pops = sorted(set(df1['pop']))

    for pop in pops:
        df = df1[df1['pop']==pop]
        controls = df[df.both_sexes==0].score
        cases = df[df.both_sexes==1].score
        print(f'{pop}: {stats.ttest_ind(controls, cases).pvalue}')
        

if __name__=='__main__':
    pops = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID']
    for logscale in [True]:
        for pop in pops:
            plot_prs_corr_hist(pop=pop, logscale=logscale)
    for pop in pops:
        plot_prs_corr_holdout(pop=pop)
    plot_score_phen(fname='100001.prs.tsv.gz', p_threshold=0.5)
    plot_score_phen(fname='100001.prs.tsv.gz', pop='AFR')
    plot_score_phen(fname='100320.prs.tsv.gz', logscale=True)
    plot_score_phen(fname='100320.prs.tsv.gz', force_quant=True)
    plot_score_phen(fname='100020.prs.tsv.gz')
    plot_score_phen(fname='100020.prs.tsv.gz', p_threshold=1)
    plot_score_phen(fname='100020.prs.tsv.gz', logscale=True)
    plot_score_phen(fname='100020.prs.tsv.gz', force_quant=True)
    plot_score_phen(fname='100020.prs.tsv.gz', pop='CSA')
    plot_score_phen(fname='250.2.prs.tsv.gz', p_threshold=1)
    plot_score_phen(fname='30600.prs.tsv.gz')
    plot_score_phen(fname='E78.prs.tsv.gz', p_threshold=1)
    
    test_prs_score_diff(fname='E78.prs.tsv.gz')
    
    plot_r2_pval_threshold(fname='100001.prs.tsv.gz', pop='AFR')
    plot_r2_pval_threshold(fname='100020.prs.tsv.gz', pop='CSA')
    plot_r2_pval_threshold(fname='250.2.prs.tsv.gz', pop='AMR')
    plot_r2_pval_threshold(fname='E78.prs.tsv.gz', pop='CSA')