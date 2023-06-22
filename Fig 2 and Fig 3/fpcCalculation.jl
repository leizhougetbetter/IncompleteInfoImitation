using Random
using Graphs, SimpleWeightedGraphs
using DelimitedFiles
using LinearAlgebra
using StatsBase




struct nodeProperty
    numNeighborsNode::Int64
    neighborsNode::Vector{Int64}
    neighborsNodeWeights::Vector{Int64}
    nodeInNeighborIndex::Vector{Int64}
    totalWeightNode::Int64
    probToSelectEachNeighborSorted::Vector{Float64}
    intevalToSelectEachNeighborSorted::Vector{Float64}
    neighborsSortedToSelectionProb::Vector{Int64}
end

function readFixedNetwork(networkFolder, networkType, weightDist)
    fileName = "NetworkStructure_" * networkType * "_Weighted_" * weightDist
    filePath =  networkFolder * '\\' * fileName
    ### read delimitedFiles with ',' as delimiter
    adjGraph = readdlm(filePath * ".dlm", ',',Int)
    #### Construct weighted Graphs
    networkImport = Graphs.SimpleGraph(adjGraph)
    ####
    NtotalNode = Graphs.nv(networkImport)
    networkStructure = Vector{nodeProperty}(undef, NtotalNode)
    if NtotalNode != Ntotal
        error("Network size does not match!")
    end
    for iNode in 1:NtotalNode
        neighborsNow = Graphs.neighbors(networkImport, iNode)
        numNeighborsNodeNow = length(neighborsNow)
        # weightsNeighborNow = networkImport.weights[neighborsNow, iNode]
        weightsNeighborNow = ones(numNeighborsNodeNow)
        nodeInNeighborIndexNow = zeros(numNeighborsNodeNow)
        totalInteractionTimesWithNeighborsNow = sum(weightsNeighborNow)
        probToSelectNeighborsNow = (weightsNeighborNow * 1.0) / totalInteractionTimesWithNeighborsNow

        ###### sort the probability to descending order (accelerating random selection)
        indexSortedToSelectionProb = sortperm(probToSelectNeighborsNow, rev=true)
        neighborsSortedToSelectionProb = neighborsNow[indexSortedToSelectionProb]
        probToSelectNeighborsSorted = probToSelectNeighborsNow[indexSortedToSelectionProb]

        intervalToSelectNeighborsNowSorted = zeros(numNeighborsNodeNow)
        cumsum!(intervalToSelectNeighborsNowSorted, probToSelectNeighborsSorted)
        pushfirst!(intervalToSelectNeighborsNowSorted, 0.0)
        ###
        networkStructure[iNode] = nodeProperty(numNeighborsNodeNow, neighborsNow, weightsNeighborNow, nodeInNeighborIndexNow, totalInteractionTimesWithNeighborsNow, probToSelectNeighborsSorted, intervalToSelectNeighborsNowSorted, neighborsSortedToSelectionProb)

    end
    #########
    for iNode in 1:NtotalNode
        degreeNodeNow = networkStructure[iNode].numNeighborsNode
        for jNode in 1: degreeNodeNow
            neighborNow = networkStructure[iNode].neighborsNode[jNode]
            neighborsOfNeighbor = copy(networkStructure[neighborNow].neighborsNode)
            networkStructure[iNode].nodeInNeighborIndex[jNode] = findfirst(xx -> xx == iNode, neighborsOfNeighbor)
        end
    end

    return networkStructure
end

function initializeStrategy(Ntotal)
    #### C  = 1, D = 0
    
    #### initialize the game state with all-D configuration
    strategyNode = zeros(Int, Ntotal);
    randNode_APlayer = rand(1:Ntotal);
    #### randomly select one individual as A-player
    strategyNode[randNode_APlayer] = 1;
    return randNode_APlayer, 1
end
#
function interactionGamesV2(nodeCalculate, strategyMain, gamePayoff, networkStructure)
    ### for pairwise games
    ### Calculate the payoff of a given node, i.e., nodeCalculate
    strategyNode = strategyMain[nodeCalculate]
    totalInteractionTimesWithNeighbors = networkStructure[nodeCalculate].totalWeightNode
    payoffNode = 0.0
    numNeighbors = networkStructure[nodeCalculate].numNeighborsNode
    payoffNode = payoffCalculationV2(payoffNode, strategyNode, numNeighbors, gamePayoff, networkStructure[nodeCalculate].neighborsNodeWeights, strategyMain, networkStructure[nodeCalculate].neighborsNode)
    ### Average payoff
    payoffNode = payoffNode / totalInteractionTimesWithNeighbors
    ### Accumulated payoff

    return payoffNode
end
#
function payoffCalculationV2(payoffNode, strategyNode, numNeighbors, gamePayoff, weightsNeighborsNode, strategyMain, neighborsNodeNow)
    ### for pairwise games
    for indNeighbor in 1 : numNeighbors
        weightNow = weightsNeighborsNode[indNeighbor]
        strategyNeighborNow = strategyMain[neighborsNodeNow[indNeighbor]]
        payoffNode += weightNow * (strategyNeighborNow * 1.0 * gamePayoff[2-strategyNode,1] + (1.0 - strategyNeighborNow) * gamePayoff[2-strategyNode,2])
    end
    return payoffNode
end
#
function interactionPGGGames(nodeCalculate, strategyMain, gamePayoff, networkStructure)
    ### for public goods games
    ### Calculate the payoff of a given node, i.e., nodeCalculate
    strategyNode = strategyMain[nodeCalculate]
    totalInteractionTimesWithNeighbors = networkStructure[nodeCalculate].totalWeightNode
    payoffNode = 0.0
    numNeighbors = networkStructure[nodeCalculate].numNeighborsNode
    ### In a PGG, gamePayoff = [bb, cc]
    ### focal PGG
    payoffNode += payoffInOnePGG(numNeighbors, gamePayoff[1], strategyNode, networkStructure[nodeCalculate].neighborsNodeWeights, strategyMain, networkStructure[nodeCalculate].neighborsNode)
    ### PGG organized by all neighbors
    for indNeighbor in 1 : numNeighbors
        neighborNow = networkStructure[nodeCalculate].neighborsNode[indNeighbor]
        strategyNeigh = strategyMain[neighborNow]
        payoffNode += payoffInOnePGG(networkStructure[neighborNow].numNeighborsNode, gamePayoff[1], strategyNeigh, networkStructure[neighborNow].neighborsNodeWeights, strategyMain, networkStructure[neighborNow].neighborsNode)
    end
    payoffNode -= (numNeighbors+1)*strategyNode*gamePayoff[2]
    ### Average payoff
    payoffNode = payoffNode / (numNeighbors+1)
    ### Accumulated payoff

    return payoffNode
end
#
function payoffInOnePGG(numNeighbors, benefitB, strategyNode, weightsNeighborsNode, strategyMain, neighborsNodeNow)
    ### for public goods games
    sumCoop = 0
    sumCoop += strategyNode
    for indNeighbor in 1 : numNeighbors
        sumCoop += strategyMain[neighborsNodeNow[indNeighbor]]
    end
    ### gamePayoff in a PGG is the [benefit b, cost c]
    benefitDistributed = 1.0/(numNeighbors+1) * sumCoop * benefitB
    return benefitDistributed
end
#
#
function sImitationUpdate(nodeChosenUpdate, numNeighChosenArray, strategyMain, gamePayoff, networkStructure, selectionIntensity, theta)
    ###
    strategyNodeChosenUpdate = strategyMain[nodeChosenUpdate]
    ####
    numNeighChosen = numNeighChosenArray[nodeChosenUpdate]
    degreeNodeChosenUpdate = networkStructure[nodeChosenUpdate].numNeighborsNode
    neighborsSetSample = sample(1:degreeNodeChosenUpdate, numNeighChosen, replace = false)

    sumCoopSelectedAndSelf = 0
    sumCoopSelectedAndSelf += strategyNodeChosenUpdate

    for iNeighbor in 1 : numNeighChosen
        neighIndexNow = neighborsSetSample[iNeighbor]
        neighborNow = networkStructure[nodeChosenUpdate].neighborsNode[neighIndexNow]
        sumCoopSelectedAndSelf += strategyMain[neighborNow]
    end

    if sumCoopSelectedAndSelf == numNeighChosen+1
        diffStrChangeDirectionNodeUpdate = 0
        return diffStrChangeDirectionNodeUpdate
    elseif sumCoopSelectedAndSelf == 0
        diffStrChangeDirectionNodeUpdate = 0
        return diffStrChangeDirectionNodeUpdate
    end


    sumCoopFitness = 0.0
    sumDefectFitness = 0.0
    
    ### for pairwise social dilemma
    # fitnessFocalPlayer = exp(selectionIntensity* interactionGamesV2(nodeChosenUpdate, strategyMain, gamePayoff, networkStructure))

    ### for group social dilemma
    fitnessFocalPlayer = exp(selectionIntensity* interactionPGGGames(nodeChosenUpdate, strategyMain, gamePayoff, networkStructure))
    
    sumCoopFitness += strategyNodeChosenUpdate * theta * fitnessFocalPlayer
    sumDefectFitness += (1-strategyNodeChosenUpdate) * theta * fitnessFocalPlayer
    
    for iNeighbor in 1 : numNeighChosen
        neighIndexNow = neighborsSetSample[iNeighbor]
        neighborNow = networkStructure[nodeChosenUpdate].neighborsNode[neighIndexNow]
        strategyNeighborNow = strategyMain[neighborNow]

        ### for pairwise social dilemma
        # fitnessNeighborNow = exp(selectionIntensity*interactionGamesV2(neighborNow, strategyMain, gamePayoff, networkStructure))
    
        ### for group social dilemma
        fitnessNeighborNow = exp(selectionIntensity*interactionPGGGames(neighborNow, strategyMain, gamePayoff, networkStructure))

        if strategyNeighborNow == 1
            sumCoopFitness += (1-theta) / numNeighChosen * fitnessNeighborNow
        elseif strategyNeighborNow == 0
            sumDefectFitness += (1-theta) / numNeighChosen * fitnessNeighborNow  
        end
    end

    totalFitnessNeighSelectedAndSefl = sumCoopFitness + sumDefectFitness

    probUpdateToC = sumCoopFitness * 1.0 / totalFitnessNeighSelectedAndSefl
    ### strategy updating if the strategies of role model and focal player differ
    randNumberNow = rand()
    if randNumberNow <= probUpdateToC
        ### imitate the strategy of the role model
        strategyMain[nodeChosenUpdate] = 1
    else
        strategyMain[nodeChosenUpdate] = 0
    end
    diffStrChangeDirectionNodeUpdate = strategyMain[nodeChosenUpdate] - strategyNodeChosenUpdate
    return diffStrChangeDirectionNodeUpdate
end
#
function calculateFPCorD(numRepetition, numGeneration, networkStructure, gamePayoff, selectionIntensity, numNeighChosen, theta)
    ####
    Ntotal = length(networkStructure)
    ###
    invasionTimes = 0
    extinctionTimes = 0
    for iRepetition in 1 : numRepetition
        if iRepetition % 5000000 == 0
            println("Repetition # = $iRepetition")
        end
        ##### Strategy and state initialization
        nodeChosenUpdate, strChangeDirectionCurrent  = initializeStrategy(Ntotal)
        strategyMain = zeros(Int, Ntotal)
        strategyMain[nodeChosenUpdate] += strChangeDirectionCurrent
        numAPlayer = 1
        ####
        for iGeneration in 1 : numGeneration
            nodeChosenUpdate = rand(1:Ntotal)
            strChangeDirectionCurrent = sImitationUpdate(nodeChosenUpdate, numNeighChosen, strategyMain, gamePayoff, networkStructure, selectionIntensity, theta)
            numAPlayer += strChangeDirectionCurrent
            if numAPlayer == Ntotal
                invasionTimes = invasionTimes + 1
                ### exit the iGeneration loop
                break
            elseif numAPlayer == 0
                extinctionTimes = extinctionTimes + 1
                break
            end
        end
    end
    fpC = (invasionTimes * 1.0) / numRepetition
    return fpC, invasionTimes, extinctionTimes, numRepetition - invasionTimes - extinctionTimes
end
#
#
function mainProgramCalculateFpCandD(pairwiseGamePars::Tuple{Float64, Float64, Float64, Float64}, invasionTrialsPars::Tuple{Int,Int,Float64}, fileIOPars::NTuple{3,String}, numNeighChosenMain, theta)
    networkFolder, networkType, weightDist = fileIOPars
    outputFileName01 = networkFolder * "\\" * "ConditionalFixationProbTime_" * networkType * "_" * weightDist
    numRepetition, numGeneration, selectionIntensityArray = invasionTrialsPars
    SSBeg, SSEnd, SSGap, TTvalue = pairwiseGamePars

    networkStructureMain = readFixedNetwork(networkFolder, networkType, weightDist)
    gamePayoff = zeros(2, 2)
    SSValue = SSBeg : SSGap : SSEnd

    for iSelectionIntensity in 1 : length(selectionIntensityArray)
        selectionIntensity = selectionIntensityArray[iSelectionIntensity]
        outputFileName = outputFileName01 * "_w" * string(selectionIntensity)
        println("outputFileName = $outputFileName")
        #### first rho_A then rho_B
        for indSS in 1 : length(SSValue)
            SSNow = SSValue[indSS]
            println("-------** SS = $SSNow, w = $selectionIntensity **-------")
            for rhoCorRhoD in 0 : 1
                ### gamePayoff = [AA, AB; BA, BB];
                if rhoCorRhoD == 0
                    ### calculate rho_A
                    #### SSNow = bb, TTvalue = cc;
                    ##for pairwise social dilemmas
                    # gamePayoff = [SSNow-TTvalue -TTvalue; SSNow 0 ];

                    ##for group social dilemmmas
                    gamePayoff = [SSNow TTvalue];

                    fpCorD, invasionTimes, extinctionTimes, coexistenceTimes = calculateFPCorD(numRepetition, numGeneration, networkStructureMain, gamePayoff, selectionIntensity, numNeighChosenMain, theta)
                    # Juno.@profiler fpCorD, invasionTimes, extinctionTimes, coexistenceTimes = calculateFPCorD(updateFunctionMain, numRepetition, numGeneration, networkStructureMain, gamePayoff, selectionIntensity, gameTransitionProb, theta)
                elseif rhoCorRhoD == 1
                    ### calculate rho_B
                    ##for pairwise social dilemmas
                    gamePayoff = [0 SSNow; -TTvalue SSNow-TTvalue];

                    ##for group social dilemmmas
                    # gamePayoff = [-SSNow -TTvalue];
                    fpCorD, invasionTimes, extinctionTimes, coexistenceTimes = calculateFPCorD(numRepetition, numGeneration, networkStructureMain, gamePayoff, selectionIntensity, numNeighChosenMain, theta)
                end
                open(outputFileName, "a") do io
                    ##### Convert large integer numbers to String form to avoid output numbers with scientific notation
                    writedlm(io, [SSNow fpCorD string(invasionTimes) string(extinctionTimes) string(coexistenceTimes)])
                end
            end
        end
    end
end
#
function assignInputArgs(inputArgsArray, networkFolder)
    inputArgsStr = readlines(networkFolder * "\\inputParameters.inp")
    inputArgsStrArray = split.(inputArgsStr, " = ")
    if length(inputArgsStrArray) != length(inputArgsArray)
        error("The number of input args is wrong")
    end
    for indInputArgs = 1 : length(inputArgsStrArray)
        inputArgsVariableNow = inputArgsStrArray[indInputArgs][1]
        inputArgsValueNow = inputArgsStrArray[indInputArgs][2]
        if indInputArgs == 1 && inputArgsVariableNow == "networkType"
            #### networkType String
            # networkType = inputArgsValueNow
            inputArgsArray[1] = string(inputArgsValueNow)
        elseif indInputArgs == 2 && inputArgsVariableNow == "weightDist"
            #### weightDist String
            # weightDist = inputArgsValueNow
            inputArgsArray[2] = string(inputArgsValueNow)
        elseif indInputArgs == 3 && inputArgsVariableNow == "numRepetitionInput"
            #### numRepetitionInput Int64
            # numRepetitionInput = parse(Int64, inputArgsValueNow)
            inputArgsArray[3] = parse(Int64, inputArgsValueNow)
        elseif indInputArgs == 4 && inputArgsVariableNow == "numGenerationInput"
            #### numGenerationInput Int64
            # numGenerationInput = parse(Int64, inputArgsValueNow)
            inputArgsArray[4] = parse(Int64, inputArgsValueNow)
        elseif indInputArgs == 5 && inputArgsVariableNow == "selectionIntensityInput"
            #### selectionIntensityInput
            # selectionIntensityInput = parse(Float64, inputArgsValueNow)
            inputArgsArray[5] = parse(Float64, inputArgsValueNow)
        elseif indInputArgs == 6 && inputArgsVariableNow == "numNeighChosenMain"
            #### numNeighChosenMain
            # numNeighChosenMain = parse(Int64, inputArgsValueNow)
            inputArgsArray[6] = parse(Int64, inputArgsValueNow)
        elseif indInputArgs == 7 && inputArgsVariableNow == "SSEstimatedCritical"
            #### SSEstimatedCritical
            # SSEstimatedCritical = parse(Float64, inputArgsValueNow)
            inputArgsArray[7] = parse(Float64, inputArgsValueNow)
        elseif indInputArgs == 8 && inputArgsVariableNow == "SSDistance"
            #### SSDistance
            # SSDistance = parse(Float64, inputArgsValueNow)
            inputArgsArray[8] = parse(Float64, inputArgsValueNow)
        elseif indInputArgs == 9 && inputArgsVariableNow == "SSGap"
            #### SSGap
            # SSGap = parse(Float64, inputArgsValueNow)
            inputArgsArray[9] = parse(Float64, inputArgsValueNow)
        elseif indInputArgs == 10 && inputArgsVariableNow == "TTvalue"
            #### TTvalue
            # TTvalue = parse(Float64, inputArgsValueNow)
            inputArgsArray[10] = parse(Float64, inputArgsValueNow)
        elseif indInputArgs == 11 && inputArgsVariableNow == "theta"
            #### theta
            # theta = parse(Float64, inputArgsValueNow)
            inputArgsArray[11] = parse(Float64, inputArgsValueNow)
        else
            error("Wrong input files!")
        end
    end
    return inputArgsArray
end

function readNumNeighorChosenEachNode(networkFolder)
    filePath =  networkFolder * "\\sNumNeighobrInput"
    ### read delimitedFiles with ',' as delimiter
    sNumNeighorNode = readdlm(filePath * ".txt", ',', Int)
    return sNumNeighorNode[:]
end

##### UpdateRule options: "pairwiseComparison", "deathBirth"
const Ntotal = 100
##### set networkFolder as path of the current working directory
networkFolder = @__DIR__
println(networkFolder)
#### Default input arg values
networkType = "RRG"
weightDist = "Homogeneous"
numRepetitionInput = 1
numGenerationInput = 1
selectionIntensityInput = 0.10
numNeighChosenMain = 1 
SSEstimatedCritical = 9.33
SSDistance = 0.00
SSGap = 0.01
TTvalue = 1.0
############# read input args
println("----------------  Input args read Begin  ----------------")
inputArgsArray = Vector(undef, 11)
inputArgsArray = assignInputArgs(inputArgsArray, networkFolder)
############# current arg values
networkType = inputArgsArray[1]
weightDist = inputArgsArray[2]
numRepetitionInput = inputArgsArray[3]
numGenerationInput = inputArgsArray[4]
selectionIntensityInput = inputArgsArray[5]
##### all indiviudals use a shared $s$
# numNeighChosenMain = inputArgsArray[6]
##### different indiviudals use different $s$
numNeighChosenMain = readNumNeighorChosenEachNode(networkFolder)
#####
SSEstimatedCritical = inputArgsArray[7]
SSDistance = inputArgsArray[8]
SSGap = inputArgsArray[9]
TTvalue = inputArgsArray[10]
theta = inputArgsArray[11]
####
SSBeg = SSEstimatedCritical - SSDistance
SSEnd = SSEstimatedCritical + SSDistance
####
println("  NetworkType = $networkType - $weightDist, networkSize = $Ntotal")
println("  numRepetitionTotal = $numRepetitionInput")
println("  numGenerationEachRep = $numGenerationInput")
println("  selectionIntensity = $selectionIntensityInput")
println("  numNeighChosen = $numNeighChosenMain")
println("  theta = $theta")
println("  TTvalue = $TTvalue (fixed)")
println("  SSDistance = $SSDistance")
println("  SSValue = $SSBeg:$SSGap:$SSEnd")
println("----------------  Input args read Done  ----------------")
#############

indGraphBeg = 1
indGraphEnd = 1
for indGraph = indGraphBeg : indGraphEnd
    println("Graph # = $indGraph")
    pairwiseGameParsInput = (SSBeg, SSEnd, SSGap, TTvalue)
    fileIOParsInput = (networkFolder, networkType, weightDist)
    ### precompiling (warm up)
    invasionTrialsParsInput = (10, 10, 1.0)
    mainProgramCalculateFpCandD(pairwiseGameParsInput, invasionTrialsParsInput, fileIOParsInput, numNeighChosenMain, theta)
    ### main function for computation
    invasionTrialsParsInput = (numRepetitionInput, numGenerationInput, selectionIntensityInput)
    mainProgramCalculateFpCandD(pairwiseGameParsInput, invasionTrialsParsInput, fileIOParsInput, numNeighChosenMain, theta)
end
println("-----------Calculation done---------------")
