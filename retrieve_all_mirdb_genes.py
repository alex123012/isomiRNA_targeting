import pandas as pd
import numpy as np
import requests as rq
from tqdm.auto import tqdm

from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi import GeneApi as DatasetsGeneApi
from ncbi.datasets.openapi.models import V1alpha1GeneMatch as V1GeneMatch

def hgnc_get_ensmbl(id_):
    headers = {
        'Accept': 'application/json',
    }

    url = f'http://rest.genenames.org/fetch/hgnc_id/{id_}'

    r = rq.get(url, headers=headers)

    return r.json()['response']['docs'][0]

def get_ensg(gene):
    gene_reply = gene_api.gene_metadata_by_accession([gene])
    ensmbl_list = set()
    synonims_list = set()

    if len((tmp := gene_reply.genes)) == 1 and not tmp[0].warnings and not tmp[0].errors:
        tmp = tmp[0].gene

        if (ensmbl := tmp.ensembl_gene_ids):
            ensmbl_list.add(ensmbl[0])
        if (synonyms := tmp.synonyms):
            for i in synonyms:
                synonims_list.add(i)
        if (symbol := tmp.symbol):
            synonims_list.add(symbol)

        if (nomen := tmp.nomenclature_authority) and nomen.authority == 'HGNC' and (hgnc := nomen.identifier) and (hgnc := hgnc_get_ensmbl(hgnc.split(':')[1])):

            if (ensmbl := hgnc.get('ensembl_gene_id')):
                ensmbl_list.add(ensmbl)
            if (synonyms := hgnc.get('prev_symbol')):
                for i in synonyms:
                    synonims_list.add(i)
            if (synonyms := hgnc.get('alias_symbol')):
                for i in synonyms:
                    synonims_list.add(i)
            if (symbol := hgnc.get('symbol')):
                synonims_list.add(symbol)

    if ensmbl_list or synonims_list:
        return {'ensmbl': list(ensmbl_list),
                'synonyms': list(synonims_list)}
    else:
        None


all_mir = pd.read_csv('/home/user503/Desktop/miRDB_v6.0_prediction_result.txt', sep='\t', header=None)
all_mir = all_mir[all_mir[0].str.split('-', 1).str[0] == 'hsa'].copy()
gene_symbols = all_mir[1].unique().tolist()

test = {}
print('not mapped symbols:', len(gene_symbols))
bad_ids = []

with DatasetsApiClient() as api_client:
    gene_api = DatasetsGeneApi(api_client)

    for gene in tqdm(gene_symbols):
        try:
            if (dictir := get_ensg(gene)):
                test[gene] = dictir

            else:
                print('no info gene:'.upper(), gene)
                bad_ids.append(gene)
                pd.DataFrame(test).T.to_csv('all_mirdb_genes.tsv')
        except Exception as exp:
            print('Error in gene:', gene)

df = pd.DataFrame(test).T
if df.ensmbl.str.len().max() == 1:
    df.ensmbl = df.ensmbl.str[0]
df.to_csv('all_mirdb_genes.tsv', sep='\t')

