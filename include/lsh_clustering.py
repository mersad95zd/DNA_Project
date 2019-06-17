'''
cluster with LSH
'''


import numpy as np
#from skbio.alignment import local_pairwise_align_ssw
#from skbio import DNA
import operator
import itertools

def kmerDNA(seq,k=4):
    kmer = []
    for ell in range(len(seq)-k+1):
        nstr = seq[ell:ell+k]
    
        index = 0
        for j,c in enumerate(nstr):
            if c == 'A':
                i = 0
            elif c == 'C':
                i = 1
            elif c == 'G':
                i = 2
            elif c == 'T':
                i = 3
            else:
                index = -1
                break
            index += i*(4**j)
        
        kmer += [ index ]
        
    return kmer


class minhashsig():
    # min hash for k-mers
    def __init__(self,m,k):
        # m is number of signatures
        self.tables = [np.random.permutation(4**k) for i in range(m)]
        self.k = k
        
    def generate_signature(self,seq):
        kmer = kmerDNA(seq,self.k)
        sig = [ min([table[i] for i in kmer])  for table in self.tables]
        return sig

    
  
        
    
#The algorithm performs a sequential scan over the sorted pairs. 
#The first time that node u appears in the scan, it is marked as a
#cluster center. All subsequent nodes v that appear in pairs
#of the form (u, v) are marked as belonging to the cluster of u 
#and are not considered again.

def center_cluster(pairs):
    clusters = {}
    for (u,v) in pairs:
        if u in clusters:
            clusters[u] += [v]
        if v in clusters:
            clusters[v] += [u]
        if v not in clusters and  u not in clusters:
            clusters[u] = [v]
    list_clusters = []
    for d in clusters:
        list_clusters += [[d] + clusters[d]]
    return list_clusters

def max_match(seq1,seq2):
    # returns number of matching characters of alignment
    alignment, score, start_end_positions \
          = local_pairwise_align_ssw(DNA( seq1 ),DNA( seq2 ), match_score=2, mismatch_score=-3)
    a = str(alignment[0])
    b = str(alignment[1])
    ctr = 0
    for i,j in zip(a,b):
        if(i==j):
            ctr +=1
    return(ctr)
#max_match("ACGTACGT","ACGTCGT")
   
    

def extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig):
    # sigs: minhash signatures
    # ell_lsh: number of LSH signatures
    # k_lsh: number of signatures to concatenate
    
    # preprocessing step
    #kmers = []
    #for seq in seqs:
    #    kmers += [ set(kmerDNA(seq,k=4)) ]
    
    pairs = set([])
    
    # generate ell_lsh random indices
    for ell in range(ell_lsh):
        # generate k_lsh independent indices from [m]
        lshinds = np.random.permutation(m)[:k_lsh]
        # generate lsh signature
        lshsigs = []
        for sig in sigs:
            lshsig = 0
            for i,lshind in enumerate(lshinds):
                lshsig += sig[lshind]*(maxsig**i)
            lshsigs += [lshsig]
        d = {}
        for ind,sig in enumerate(lshsigs):
            if sig in d:
                d[sig] += [ind]
            else:
                d[sig] = [ind]
        for candidates in d.values():
            if len(candidates) > 1:
                for pair in itertools.combinations(candidates,2):
                    # filter here
                    #score = len(kmers[u] & kmers[v]) / min(len(kmers[u]), len(kmers[v])  )
                    #if(score >= maxmatch):
                    pairs.add(pair)
    return pairs         

def filter_pairs(pairs,seqs,threshold=0.7):
    kmers = []
    for seq in seqs:
        kmers += [ set(kmerDNA(seq,k=4)) ]
    
    filtered_pairs = []
    for (u,v) in pairs:
        score = len(kmers[u] & kmers[v]) / min(len(kmers[u]), len(kmers[v])  )
        if(score >= threshold):
            filtered_pairs.append( (u,v) )
    return filtered_pairs

def lsh_cluster(seqs,m,k,k_lsh=3,ell_lsh=4,threshold=0.7):
    '''
    Unrelated strings are unlikely to agree on all  kk  min-hash signatures, 
    thus larger  kk  reduces the number of false positives. On the other hand,
    it also increases the number of false negatives. To reduce the latter effect,
    ℓℓ  different LSH signatures are extracted for each URL. This increases 
    the chance that at least two related URLs agree on at least one of their LSH signature.
    '''
    print("generate signatures")
    minhash = minhashsig(m,k)
    sigs = [minhash.generate_signature(seq) for seq in seqs]
    #print(sigs)
    maxsig = 4**k
    #print(maxsig)
    print("extract similar pairs")
    pairs = extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig)
    print("\t number pairs:",len(pairs))
    if maxmatch > 0:
        print("filter pairs")
        pairs = filter_pairs(pairs,seqs,threshold)
    print("center cluster")   
    return center_cluster(pairs)    

'''

def extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig):
    # sigs: minhash signatures
    # ell_lsh: number of LSH signatures
    # k_lsh: number of signatures to concatenate
    
    pairs = set([])
    
    # generate ell_lsh random indices
    for ell in range(ell_lsh):
        # generate k_lsh independent indices from [m]
        lshinds = np.random.permutation(m)[:k_lsh]
        # generate lsh signature
        lshsigs = []
        for sig in sigs:
            lshsig = 0
            for i,lshind in enumerate(lshinds):
                lshsig += sig[lshind]*(maxsig**i)
            lshsigs += [lshsig]
        d = {}
        for ind,sig in enumerate(lshsigs):
            if sig in d:
                d[sig] += [ind]
            else:
                d[sig] = [ind]
        for candidates in d.values():
            if len(candidates) > 1:
                for pair in itertools.combinations(candidates,2):
                    pairs.add(pair)
    return pairs    
    
def filter_pairs(pairs,seqs,maxmatch=40):
    filtered_pairs = []
    for (u,v) in pairs:
        score = max_match( seqs[u], seqs[v] )
        if(score >= maxmatch):
            filtered_pairs  +=  [(u,v)]
    return filtered_pairs

def lsh_cluster(seqs,m,k,k_lsh=2,ell_lsh=4,maxmatch=0):
    print("generate signatures")
    minhash = minhashsig(m,k)
    sigs = [minhash.generate_signature(seq) for seq in seqs]
    #print(sigs)
    maxsig = 4**k
    #print(maxsig)
    print("extract similar pairs")
    pairs = extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig)
    print("\t number pairs")
    if maxmatch > 0:
        pairs = filter_pairs(pairs,seqs,maxmatch)
    print("center cluster")   
    return center_cluster(pairs)
'''





