#!/usr/bin/env python

"""
Currently, there are issues when trying to run in Python 3, so we need to run it
in Python 2 for now
"""
import numpy as np
import pandas as pd
from scipy import stats
from numpy import random as rd
import os
import json
from sklearn.decomposition import PCA
import multiprocessing, multiprocessing.pool
import matplotlib.pyplot as plt
import time
from collections import Counter
from sklearn.metrics import f1_score
from miner import miner

RESULTS_DIR = "results"

if __name__ == '__main__':
    if not os.path.isdir(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)

    # read gene expression into pandas dataframe
    exp_data = pd.read_csv(os.path.join("minerData", "IA12Zscore.csv"),
                           index_col=0, header=0)
    #exp_data = pd.read_csv(os.path.join("..", "data", "expressionData.csv"),
    #                       index_col=0, header=0)

    # convert gene names to Ensembl Gene ID
    exp_data = miner.identifierConversion(exp_data)

    # pre-process expressionData for optimized performance
    exp_data = miner.zscore(exp_data)

    # Step 1. cluster expression data
    # generate a list of coexpressed gene clusters, all of which have length >= minNumberGenes
    init_clusters = miner.cluster(exp_data, minNumberGenes=8, overExpressionThreshold=80)

    # sort cluster list from largest to smallest
    init_clusters.sort(key= lambda s: -len(s))

    # revise initialClusters to combine highly correlated clusters and keep only
    # those with significant coexpression
    revised_clusters = miner.reviseInitialClusters(init_clusters, exp_data)
    print(revised_clusters)

    # write revisedClusters to .json file
    miner.write_json(revised_clusters,
                     os.path.join(RESULTS_DIR, "coexpressionDictionary.json"))


    # Infer bicluster membership of revised clusters
    # create a background matrix used for statistical hypothesis testing
    bkgd = miner.backgroundDf(exp_data)

    # for each cluster, give samples that show high coherent cluster activity
    overexpressed_members = miner.biclusterMembershipDictionary(revised_clusters,
                                                                bkgd, label=2)

    # for each cluster, give samples that show low coherent cluster activity
    underexpressed_members = miner.biclusterMembershipDictionary(revised_clusters,
                                                                 bkgd, label=0)

    # for each cluster, give samples that do not show coherent cluster activity
    dysregulated_members = miner.biclusterMembershipDictionary(revised_clusters,
                                                               bkgd, label="excluded")

    # for each cluster, give samples that show coherent cluster activity,
    # despite magnitude of expression
    coherent_members = miner.biclusterMembershipDictionary(revised_clusters,
                                                           bkgd, label="included")


    # Generate bicluster membership incidence matrices
    # write membership matrices to .csv file
    overexpressed_members_matrix = miner.membershipToIncidence(overexpressed_members,
                                                               exp_data)
    overexpressed_members_matrix.to_csv(os.path.join(RESULTS_DIR, "overExpressedMembers0.csv"))

    underexpressed_members_matrix = miner.membershipToIncidence(underexpressed_members,
                                                                exp_data)
    underexpressed_members_matrix.to_csv(os.path.join(RESULTS_DIR, "underExpressedMembers0.csv"))

    dysregulated_members_matrix = miner.membershipToIncidence(dysregulated_members,
                                                              exp_data)
    dysregulated_members_matrix.to_csv(os.path.join(RESULTS_DIR, "dysregulatedMembers0.csv"))

    coherent_members_matrix = miner.membershipToIncidence(coherent_members,
                                                          exp_data)
    coherent_members_matrix.to_csv(os.path.join(RESULTS_DIR, "coherentMembers0.csv"))

    cluster_median_matrix = miner.coexpressionMedianMatrix(exp_data,
                                                           revised_clusters)
    #print(cluster_median_matrix.head(5))

    # 2. Infer coregulated modules

    # get first principal component axes of clusters
    axes = miner.principalDf(revised_clusters, exp_data, subkey=None, minNumberGenes=1)

    # analyze revised clusters for enrichment in relational database
    # (default: transcription factor binding site database)
    mechanistic_output = miner.mechanisticInference(axes,
                                                    revised_clusters,
                                                    exp_data,
                                                    dataFolder=os.path.join("data"))

    # write mechanistic output to .json file
    miner.write_json(mechanistic_output, os.path.join(RESULTS_DIR, "mechanisticOutput.json"))

    # order mechanisticOutput as {tf:{coexpressionModule:genes}}
    coregulation_modules = miner.getCoregulationModules(mechanistic_output)

    # write coregulation modules to .json file
    miner.write_json(coregulation_modules,
                     os.path.join(RESULTS_DIR, "coregulationModules.json"))

    # get final regulons by keeping genes that requently appear coexpressed
    # and associated to a common regulator
    regulons = miner.getRegulons(coregulation_modules,
                                 minNumberGenes=5,
                                 freqThreshold = 0.333)

    # write regulons to json file
    miner.write_json(regulons,
                     os.path.join(RESULTS_DIR, "regulons.json"))

    # define coexpression modules as composite of coexpressed regulons
    coexpression_modules = miner.getCoexpressionModules(mechanistic_output)

    # write coexpression modules to .json file
    miner.write_json(coexpression_modules,
                     os.path.join(RESULTS_DIR, "coexpressionModules.json"))

    # Infer bicluster membership of coexpression modules
    # create a background matrix used for statistical hypothesis testing
    bkgd = miner.backgroundDf(exp_data)

    # for each cluster, give samples that show high coherent cluster activity
    overexpressed_members = miner.biclusterMembershipDictionary(coexpression_modules,
                                                                bkgd, label=2)

    # for each cluster, give samples that show low coherent cluster activity
    underexpressed_members = miner.biclusterMembershipDictionary(coexpression_modules,
                                                                 bkgd, label=0)

    # for each cluster, give samples that do not show coherent cluster activity
    dysregulated_members = miner.biclusterMembershipDictionary(coexpression_modules,
                                                               bkgd, label="excluded")

    # for each cluster, give samples that show coherent cluster activity,
    # despite magnitude of expression
    coherent_members = miner.biclusterMembershipDictionary(coexpression_modules,
                                                           bkgd, label="included")


    # Generate bicluster membership incidence matrices
    # write membership matrices to .csv file
    overexpressed_members_matrix = miner.membershipToIncidence(overexpressed_members,
                                                               exp_data)
    overexpressed_members_matrix.to_csv(os.path.join(RESULTS_DIR,"overExpressedMembers.csv"))

    underexpressed_members_matrix = miner.membershipToIncidence(underexpressed_members,
                                                                exp_data)
    underexpressed_members_matrix.to_csv(os.path.join(RESULTS_DIR, "underExpressedMembers.csv"))

    dysregulated_members_matrix = miner.membershipToIncidence(dysregulated_members,
                                                              exp_data)
    dysregulated_members_matrix.to_csv(os.path.join(RESULTS_DIR, "dysregulatedMembers.csv"))

    coherent_members_matrix = miner.membershipToIncidence(coherent_members, exp_data)
    coherent_members_matrix.to_csv(os.path.join(RESULTS_DIR, "coherentMembers.csv"))

    coexpression_matrix = miner.coexpressionMedianMatrix(exp_data, coexpression_modules)

    # -----------------------------------------------------------------------
    # Predict sample classes
    sample_dictionary = overexpressed_members
    sample_matrix = overexpressed_members_matrix

    # "fequencyThreshold" is the actual parameter name of the function
    similarity_clusters = miner.getSimilarityClusters(sample_dictionary,
                                                      fequencyThreshold=0.333,
                                                      similarityThreshold=0.15,
                                                      highResolution=False)

    hires_clusters, unclustered = miner.getSimilarityClusters(sample_dictionary,
                                                              fequencyThreshold=0.333,
                                                              similarityThreshold=0.15,
                                                              highResolution=True)

    initial_classes = miner.classesFromClusters(hires_clusters)
    sample_freq_matrix = miner.sampleCoincidenceMatrix(sample_dictionary,
                                                       freqThreshold=0.333,
                                                       frequencies=True)

    similarity_matrix = sample_freq_matrix * sample_freq_matrix.T

    plt.imshow(similarity_matrix.loc[np.hstack(initial_classes),
                                     np.hstack(initial_classes)],
               vmax=0.5)
    plt.savefig(os.path.join(RESULTS_DIR,"patientClasses15percent.pdf"))

    centroid_clusters, centroid_matrix = miner.centroids(initial_classes,
                                                         sample_matrix,
                                                         f1Threshold=0.3,
                                                         returnCentroids=True)

    centroid_matrix.to_csv(os.path.join(RESULTS_DIR, "centroids.csv"))
    unmapped = list(sample_matrix.columns[np.where(sample_matrix.sum(axis=0)==0)[0]])
    mapped_samples = [i for i in np.hstack(centroid_clusters) if i not in unmapped]
    plt.imshow(similarity_matrix.loc[mapped_samples, mapped_samples],
               vmin=0, vmax=0.55)
    plt.savefig(os.path.join(RESULTS_DIR,"mapped_samples.pdf"))


    mapped_clusters = miner.mapExpressionToNetwork(centroid_matrix,
                                                   sample_matrix,
                                                   threshold=0.25)

    ordered_overexpressed_members = miner.orderMembership(centroid_matrix,
                                                          sample_matrix,
                                                          mapped_clusters,
                                                          resultsDirectory=RESULTS_DIR)


    startTimer = time.time()
    laplacian, eigenvectors, eigenvalues = miner.laplacian(regulons)
    stopTimer = time.time()
    print('completed topological analysis in {:.2f} minutes'.format((stopTimer-startTimer)/60.))
