import sys,csv, pandas
#Library just for pretty print
from pandas import *
def alignSequences(sequenceFile):
 
    #loop through the rows of the sequence file
    for csv_row in range(sequenceFile.shape[0]):
        sequence1 = sequenceFile["sequence1"][csv_row]
        sequence2 = sequenceFile["sequence2"][csv_row]
        matrix = [[0]]
        d = -2
        #Generate 2d array with rows represented as sequence1, and columns as sequence2.  
        def initMatrix():
            for columnIndex in range(len(sequence1)):
                #Create first row
                matrix[0].append(d*(columnIndex+1))
            for rowIndex in range(len(sequence2)):
                #create first column
                matrix.append([d*(rowIndex+1)])

        def scoringMatrix(row, column):
            # if len(sequence1)-1 < row or len(sequence2)-1 < column:
            #     return None
            #top left diagonal plus the scoring matrix
            tld = matrix[row-1][column-1]
            sm = 0
            if sequence2[row-1] == sequence1[column-1]:
                sm = 1
            else:
                sm = -1
            return tld + sm
            
        def maxScore(row, column):
            #return the largest value out of the three criteria
            scoring_matrix = scoringMatrix(row, column)
            #gap penalty
            #compare with left
            gp_1 = matrix[row][column-1] + d
            #compare with top
            gp_2 = matrix[row-1][column] + d 

            #return the greatest number from the three
            # if scoring_matrix == None:
            #      return max([gp_1, gp_2])
            #else
            return max([scoring_matrix, gp_1, gp_2])
           

        initMatrix()
        for rowIndex in range(1, len(sequence2)+1):
            for columnIndex in range(1, len(sequence1)+1): 
                matrix[rowIndex].append(maxScore(rowIndex,columnIndex))
        #Proteins are the individual symbols in a sequence
        #Loop through the proteins in the first sequence
        
        #Backtracking
        #Check if sequence row and column index are equal, find gaps
        #Start at last row, last column and backtrack
        alignmentRow = ""
        alignmentColumn = ""
        row = len(sequence2)
        column = len(sequence1)
        #Keep running until bound is reached
        while row != 0 and column != 0:
            #if sequence proteins match, go up diagonally
            proteinRow = sequence2[row-1]
            proteinColumn = sequence1[column-1]

            #check if match, mismatch, or gap
            #Scores
            diag = scoringMatrix(row,column)
            # add gap penalties
            upper = matrix[row-1][column] + d
            side = matrix[row][column-1] + d

            if (diag  > upper and diag  > side):

                #Add protein to the sequence (protein1 and protein2 are the same)
                alignmentRow = proteinRow + alignmentRow 
                alignmentColumn = proteinColumn + alignmentColumn
                row = row-1
                column = column -1

            elif upper > side:
                #top side is greatest, move up and add gap 
                alignmentRow = proteinRow + alignmentRow 
                alignmentColumn = "-" + alignmentColumn

                row = row-1   
            else:
                #left side is greatest, move left and add gap
                alignmentRow = "-" + alignmentRow 
                alignmentColumn = proteinColumn + alignmentColumn
                column = column-1

        print(alignmentRow)  
        print(alignmentColumn)
      
        #pretty print matrix:
        print(DataFrame(matrix))
        #[0, A, B, C , D]
        #[A             ]
        #[B             ]
        #[C             ]
        #[D             ]

def main():
    #collect file directory as argument
    filename = sys.argv[1]
    #read the file if it exists
    csvFile = open(filename, "r")
    #Separate each sequence by the respective headers
    col_list = ["sequence1", "sequence2"]
    #create the sequence file, organized data
    sequenceFile = pandas.read_csv(csvFile, usecols=col_list)
    #print(sequenceFile["sequence2"])
    alignSequences(sequenceFile)

if __name__ == "__main__":
    main()