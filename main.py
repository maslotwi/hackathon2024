from time import sleep


def KMP_T(W):
    T = [0 for _ in range(len(W))]
    i = 2
    j = 0
    T[0] = -1
    T[1] = 0
    while i < len(W):
        if W[i - 1] == W[j]:
            T[i] = j + 1
            i = i + 1
            j = j + 1
        else:
            if j > 0:
                j = T[j]
            else:
                T[i] = 0
                i = i + 1
    return T


def KMP_search(S, W, tolerance):
    m = 0
    i = 0
    T = KMP_T(W)
    mismatch = 0
    saved = 0
    while m + i < len(S):
        if W[i] == S[m + i]:
            i += 1
            if i == len(W):
                return m, mismatch
        else:
            if mismatch == 0:
                saved = i
            mismatch += 1
            if mismatch > tolerance:
                m = m + saved - T[saved] + 1
                if i > 0:
                    i = T[saved] - 1
                mismatch = 0
            i += 1
            if i == len(W):
                return m, mismatch
    return None


# print('\n'+W[:i+1], S[m:m+i+1], mismatch_, saved, i+1)
# print(W, S)
# print(" "*i+"^"+" "*(len(W)-i)+" "*(m+i)+"^")
# print()
# if m%10000 == 0:
#     print('\r'+str(m),end="")
# sleep(1)


#
#
# text = b"GACTCAAGAGCCACGGGTGACCAGCGGCTACCGTCATGGA"
# pattern = b".AAGAG"


text = open("challenge_data/22q12_fragment/Patient_1_genome_fragment.fa", 'rb').read().replace(b'\n', b'')
pattern = open('challenge_data/anemia_gene/gene.fna', 'rb').read().replace(b'\n', b'')

result, mismatch_ = KMP_search(
    text,
    pattern,
    len(pattern) // 10
)

print()
print(result, f'{mismatch_}/{len(pattern)} ({(1-mismatch_ / len(pattern)) * 100}% match)')

# print(KMP_T(open('challenge_data/anemia_gene/gene.fna').read()))
