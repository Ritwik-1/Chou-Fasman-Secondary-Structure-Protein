# Question 1. Write a program to implement Chou and Fasman method of secondary structure
# prediction (Helix and Beta strand) for the following protein sequence.
# SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDL
# LIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLA
# CYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF
# a) Provide a well commented code for the same. (30 marks)
# b) Display the output using the notation – H: Helix, S: Beta strand. (10 marks)
# Chou and Fasman parameters to be used for the prediction are as following:


# Algorithm for Helix:-
# Take two pointers i and j , j=i+6 and then i++ till 
# we fins a window with propensity scores are >= 4 or j >= len(string)
# basically see score(string[i:j] >= 4)
# If we find a window with score >=4 then we check extensions both on the 
# left and right i.e. j+1,j,j-1,j-2 amino acids must have score >=4 for i+1 extension
# and i-1,i,i+1,i+2 must have score >= 4 for i-1 extension

# Algorithm for Beta sheet:
# same as helix , here we take a window of 5 and select 
# a window with score >= 3 for nucleation
# then for extension we take score >=4
# propensities are taken from the beta sheet dictionary


def calculateScoreHelix(s):
    score = 0
    for i in s:
        score += propensities_helix[i]
    
    return score

def calculateScoreBeta(s):
    score = 0
    for i in s:
        score += propensities_beta[i]
    
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
    if(calculateScoreHelix(sequence[i:j]) >= 4):
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
    if(calculateScoreBeta(sequence[i:j]) >= 3):
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
        
            
print("Answer with helices : ","".join(answer_helix))
print("Answer with beta sheets : ","".join(answer_beta))


print("Now resolving common regions :-")

# This will store the final answer, i.e final secondary structure 
# after accounting the conflicting regions

final_ans = list(sequence)


# For conflicting regions we will take the score of the
# sequence of conflict and see where it is greater ,i.e 
# whether in case of B sheets or alpha helix

# For this we will start i=0 and continue i++
# till any one sequence is not H or B
# then for i to i+k , k being length of common region
# we will calculate scores and assign accordingly 
  

i = 0


while i < len(answer_helix):
      if(answer_helix[i] == "H" and answer_beta[i] == "B"):
          # If the amino acids are H and B on both 
          # Then check till where are they same
          start = i
          while(answer_helix[i] == "H" and answer_beta[i] == "B"):
              i+=1
          end = i

          # when we get the start and end of the conflicting region
          # then calculate score
          score_beta = calculateScoreBeta(sequence[start:end])
          score_helix = calculateScoreHelix(sequence[start:end])

          print("Common region starts at index : ",start)
          print("Common region ends at index : ",end)
          print("score_beta : ",score_beta)
          print("score_helix : ",score_helix)
            
          # assign the conflicting region based on score answers
          if(score_beta > score_helix):
              final_ans[start:end] = "B"*((end-start))
          else:
              final_ans[start:end] = "H"*((end-start))
      else:
          # If the Beta sequence is not B and Helix is H
          # Then the secondary structure must be H at that position
          if(answer_helix[i] == "H" and answer_beta[i] != "B"):
              final_ans[i] = "H"
              i+=1
          # If the Beta sequence is B and Helix is not H 
          # Then the secondary structure must be B at that position
          elif(answer_beta[i] == "B" and answer_helix[i] != "H"):
              final_ans[i] = "B"
              i+=1
          # If none of the secondary structure can be assigned then we 
          # add the original sequence as it is
          else:
              final_ans[i] = sequence[i]
              i+=1

print("Final answer : ","".join(final_ans))
