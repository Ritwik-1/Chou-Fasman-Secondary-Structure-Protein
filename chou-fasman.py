# Question 1. Write a program to implement Chou and Fasman method of secondary structure
# prediction (Helix and Beta strand) for the following protein sequence.
# SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDL
# LIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLA
# CYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF
# a) Provide a well commented code for the same. (30 marks)
# b) Display the output using the notation – H: Helix, S: Beta strand. (10 marks)
# Chou and Fasman parameters to be used for the prediction are as following:


# For H : take two pointers i and j , j=i+6 and then i++ j++ till propensity scores are >= 4  
# or j >= len(string)
# basically see score(string[i:j] >= 4)
#if that is true then store this score and check for the element i-- and j++ if true then add 
#else i = j+1 and j = i+6

def calculateScore(s):
    score = 0
    for i in s:
        score += propensities_helix[i]
    
    return score

def checkIfcanAddToright(k):
    ele1 = propensities_helix[sequence[k]]
    ele2 = propensities_helix[sequence[k-1]]
    ele3 = propensities_helix[sequence[k-2]]
    ele4 = propensities_helix[sequence[k-3]]

    if(ele1 + ele2 + ele3 + ele4 >= 4):
        return 1
    return 0

def checkIfcanAddToLeft(k):
    # taking k = i-1
    ele1 = propensities_helix[sequence[k]]
    ele2 = propensities_helix[sequence[k+1]]
    ele3 = propensities_helix[sequence[k+2]]
    ele4 = propensities_helix[sequence[k+3]]

    if(ele1 + ele2 + ele3 + ele4 >= 4):
        return 1
    return 0

propensities_helix = {"E":1.53,
                      "A":1.45,
                      "L":1.34,
                      "H":1.24,
                      "M":1.20,
                      "Q":1.17,
                      "W":1.14,
                      "V":1.14,
                      "F":1.12,
                      "K":1.07,
                      "I":1.00,
                      "D":0.98,
                      "T":0.82,
                      "S":0.79,
                      "R":0.79,
                      "C":0.77,
                      "N":0.73,
                      "Y":0.61,
                      "P":0.59,
                      "G":0.53}

propensities_beta = { "M":1.67,
                      "V":1.65,
                      "I":1.60,
                      "C":1.30,
                      "Y":1.29,
                      "F":1.28,
                      "Q":1.23,
                      "L":1.22,
                      "T":1.20,
                      "W":1.19,
                      "A":0.97,
                      "R":0.90,
                      "G":0.81,
                      "D":0.80,
                      "K":0.74,
                      "S":0.72,
                      "H":0.71,
                      "N":0.65,
                      "P":0.62,
                      "E":0.26}

sequence = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"


# these are the answers which initially store the sequence only and
# The code will replace the values to H or B
answer_helix = sequence
answer_beta = sequence

answer_helix = list(answer_helix)
answer_beta = list(answer_beta)


i = 0
j = i+6

while(j < len(sequence)):
    if(calculateScore(sequence[i:j]) >= 4):
        for e in range(i,j):
            answer_helix[e] = "H"
        # now check for extensions on both sides 
        # right side
        k1 = j
        k2 = i-1
        #k1 is the pointer to check right extension
        while(checkIfcanAddToright(k1) == 1 and k1<len(sequence)):
            answer_helix[k1] = "H"
            k1+=1
        # now checking if can extend to the left 
        while(checkIfcanAddToLeft(k2) == 1 and k2>0):
            answer_helix[k2] = "H"
            k2-=1

        # here code will reach when k2 and k1 cannot be added to Helices
        # so will do i = k1 and j = i+6 and repeat the outermost while loop

        i = k1
        j = i + 6
    else:
        # if the score is not > 4 then we shift the window of size 6 to 
        # 1 position to the right
        i+=1
        j+=1


# For beta sheet 
i = 0
j = i+5
while(j < len(sequence)):
    if(calculateScore(sequence[i:j]) >= 3):
        for e in range(i,j):
            answer_beta[e] = "B"
        # now check for extensions on both sides 
        # right side
        k1 = j
        k2 = i-1
        #k1 is the pointer to check right extension
        while(checkIfcanAddToright(k1) == 1 and k1<len(sequence)):
            answer_beta[k1] = "B"
            k1+=1
        # now checking if can extend to the left 
        while(checkIfcanAddToLeft(k2) == 1 and k2>0):
            answer_beta[k2] = "B"
            k2-=1

        # here code will reach when k2 and k1 cannot be added to Helices
        # so will do i = k1 and j = i+6 and repeat the outermost while loop

        i = k1
        j = i + 5
    else:
        # if the score is not > 4 then we shift the window of size 6 to 
        # 1 position to the right
        i+=1
        j+=1
        
            
print("".join(answer_helix))
print("".join(answer_beta))



# Doubts : Do we check for left side extension as well in beta sheets
#          Is the nucleation site in beta >=3
#          and for extension is it >=4 or >4 
# finally we check for that the common region or at each window?
# if i -> i+6 we get helix do we check for beta sheet there only ? and change 