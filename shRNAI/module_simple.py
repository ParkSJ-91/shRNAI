import numpy as np

def pair(seq1,seq2,index):
    # make pair information between upper and lower strand of pri-miRNA
    # includes U-G pair
    if index < 12:
        return 0
    else:
        if seq1 == "A" and seq2 == "T":
            return 1
        elif seq1 == "T" and seq2 in ["A","G"]:
            return 1
        elif seq1 == "C" and seq2 == "G":
            return 1
        elif seq1 == "G" and seq2 in ["C","T"]:
            return 1
        else:
            return 0

def convert(seq):
    # make inputs from sequence
    # both mature-shRNA and pri-shRNA
    trans = str.maketrans('ACGT','0123')
    trans2 = str.maketrans('ACGT','TGCA')

    seqK = []; priK = []; onehotK = []; onehotK_pri = []
    for i in range(len(seq) - 21):
        pRNA = seq[i:i+22]
        gRNA = pRNA.translate(trans2)[::-1]
        # gRNA's 5'end should be unpaired
        if gRNA[-1] in ['A','T']:
            pRNA = 'C' + pRNA[1:]
        elif gRNA[-1] in ['C','G']:
            pRNA = 'A' + pRNA[1:]
        if gRNA in seqK:
            continue
        seqK.append(gRNA)

        var = list(map(int, list(gRNA.translate(trans))))
        var_onehot = np.eye(4)[var]
        onehotK.append(var_onehot)

        # pri-miR-30 templates
        seq_pri = "GGTATATTGCTGTTGACAGTGAGCG" + pRNA + "TAGTGAAGCCACAGATGTA" + gRNA + "TGCCTACTGCCTCGGAATTCAAGGG"
        priK.append(seq_pri)
        var_pri = list(map(int, list(seq_pri.translate(trans))))
        tempIn1 = np.delete(np.eye(5)[var_pri[:56]],-1,1)
        tempIn2 = np.delete(np.eye(5)[var_pri[-56:][::-1]],-1,1)
        tempIn3 = []
        for i in range(len(tempIn1)):
            seq1 = seq_pri[:56][i] #tempIn1[i]
            seq2 = seq_pri[-56:][::-1][i]#tempIn2[i]
            tempIn3.append(pair(seq1,seq2,i))

        seqIn = np.append(np.append(tempIn1, tempIn2, axis=1), np.asarray(tempIn3).reshape(len(tempIn3),1), axis=1)
        onehotK_pri.append(seqIn)

    onehotK = np.asarray(onehotK).reshape(-1,22,4,1)
    onehotK_pri = np.asarray(onehotK_pri).reshape(-1,56,9,1)
    return seqK, priK, onehotK, onehotK_pri

def get_Annotation(annoF):
    # grep gene annotations 
    # gencode v36 annotation
    annoDic = dict()
    f = open(annoF); lines = f.readlines(); f.close()
    for line in lines:
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        line_type = line[2]
        if line_type != 'transcript': continue

        info_line = line[8].strip(';').split('; ')
        symbol = list(filter(lambda x: x.split(' ')[0] == 'gene_name', info_line))[0].split('"')[1]
        txnID = list(filter(lambda x: x.split(' ')[0] == 'transcript_id', info_line))[0].split('"')[1]
        if txnID in ['ENST00000615113.1']: continue
        if not symbol in annoDic: annoDic[symbol] = []
        annoDic[symbol].append(txnID)
    return annoDic

def get_Sequence(inF, region, annoDic):
    # grep gene sequence
    # gencode v36 annotation
    seqDic = dict(); seqDic2 = dict(); pairDic = dict()
    f = open(inF); all = f.read(); f.close()
    chunks = all.split('>')
    for chunk in chunks[1:]:
        lines = chunk.split('\n')

        info_line = lines[0].split('|')

        txnID = info_line[0]
        geneID = info_line[1]
        symbol = info_line[5]
        if not txnID in annoDic[symbol]: continue
        cdsCoord = list(filter(lambda x: x.split(':')[0] == 'CDS', info_line))
        if len(cdsCoord) != 1:
            print(info_line)
            continue
        cdsStart = int(cdsCoord[0].split(':')[1].split('-')[0])
        cdsEnd = int(cdsCoord[0].split(':')[1].split('-')[1])

        seq = ''.join(lines[1:])
        if region == 'CDS':
            seqDic[txnID] = seq[cdsStart-1:cdsEnd]
        elif region == 'UTR':
            seqDic[txnID] = seq[cdsEnd:]
        if not symbol in pairDic: pairDic[symbol] = []
        pairDic[symbol].append(txnID)
    return seqDic, pairDic

