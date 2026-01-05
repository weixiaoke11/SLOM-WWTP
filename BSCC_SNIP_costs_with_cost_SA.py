# coding=GBK
# ======================================================================================
# Copyright 2014  Swiss Federal Institute of Aquatic Science and Technology
#
# This file is part of SNIP (Sustainable Network Infrastructure Planning)
# SNIP is used for determining the optimal degree of centralization for waste
# water infrastructures. You find detailed information about SNIP in Eggimann et al. (2014).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
#    
# The Djikstra and a* algorithm are adapted from Hetland (2010).
# The algorithm is developed for Python 2.7 and ArcGIS 10.2
#
# Literature
# ----------
# Eggimann Sven, Truffer Bernhard, Maurer Max (2015): To connect or not to connect? 
# Modelling the optimal degree of centralisation for wastewater infrastructures.   
# Water Research, XY, .....
#
# Hetland M.L. (2010): Python Algorithms. Mastering Basic Algorithms in the Python Language. apress.
#
# Contact:   sven.eggimann@eawag.ch
# Version    1.0
# Date:      1.1.2015
# Autor:     Eggimann Sven
# ======================================================================================
import math


def calculatePipeCosts(distance, Q, pipeDiameter, slope, fc_SewerCost, coefficients):
    """
    this function calculates the CH4-equal from pipelines:
        GWP_CH4 = 56 (100years)
        CH4PerDayPerKm = 0.419 * 1.06**(T-20) * Q**0.26 * D**0.28 * S**-0.138
        in which: CH4PerDayPerKm: [kgCH4/(km*day)]
                  T = 20[degree Celsius]
                  Q : flow[m3/s]
                  D : pipe diameter [m]
                  S : pipe slope [m/m]
        therefore unit conversion:
            L = distance/1000
            Q = Q/86400
            D = pipeDiameter
            S = slope
        CO2PerYear = GWP_CH4 * 365[day] * L * CH4PerDayPerKm
                   = 56 * 365 * (distance/1000) * 0.419 * 1 * (Q/86400)**0.26 * pipeDiameter**0.28 * slope**-0.138
                   = 56 * 365 * 0.419 * 86400**-0.26 / 1000 * distance * Q**0.26 * pipeDiameter**0.28 * slope**-0.138
                   = 0.44586 * distance * Q**0.26 * pipeDiameter**0.28 * slope**-0.138
        if GWP_CH4=26(20years), 0.44586 -> 0.207

    :param distance: Pipe Distance [m]
    :param Q: Flow [m3/d]
    :param pipeDiameter: m
    :param slope: %
    :param fc_SewerCost: sensitive analyze
    :return:
    GhgSewerPerYear [kgCO2/year]
    CostSewerPerYear [10^4 CNY/year]
    """
    sensFactor_b = 1 + fc_SewerCost
    CH4SewerPerYear = 0.44586 * distance * (Q ** 0.26) * (pipeDiameter ** 0.28) * (abs(slope) ** -0.138) * sensFactor_b
    if pipeDiameter < 0.6:
        CostSewerPerYear = (coefficients["Pipe1"] * distance/1000)/10000
        GhgSewerPerYear = coefficients["Pipe1_ghg"] * distance/1000 + CH4SewerPerYear
    elif pipeDiameter < 1:
        CostSewerPerYear = (coefficients["Pipe2"] * distance/1000)/10000
        GhgSewerPerYear = coefficients["Pipe2_ghg"] * distance/1000 + CH4SewerPerYear
    elif pipeDiameter < 1.5:
        CostSewerPerYear = (coefficients["Pipe3"] * distance/1000)/10000
        GhgSewerPerYear = coefficients["Pipe3_ghg"] * distance/1000 + CH4SewerPerYear
    else:
        CostSewerPerYear = (coefficients["Pipe4"] * distance/1000)/10000
        GhgSewerPerYear = coefficients["Pipe4_ghg"] * distance/1000 + CH4SewerPerYear

    return GhgSewerPerYear, CostSewerPerYear,


def getPumpCostsDependingOnFlow(Q, heightDifference, coefficients):
    """
    This function calculates the costs of a individual pump:
        gravity = 9.81[m / s^2]
        DaysPerYear = 365
        efficiency = 0.7
        motorPowerInput = (gravity * Q * heightDifference)/(efficiency)
                        = (9.81[N/kg] * Q[10^3 kg/d] * heightDifference[m])/0.7
                        = 9.81*Q*heightDifference / 0.7 kJ/d
                        = 9.81*Q*heightDifference / (0.7*3600) kWh/d
                        = 0.00389285714286 * Q*heightDifference
        EnergyUsed = motorPowerInput * DaysPerYear

    :param Q: Flow [m3/d]=[10^3 kg/d]
    :param heightDifference: Slope
    :param GridEmiFactor: grid emission factor
    :return:
    PumpGhgPerYear  GHG emissions from pumps per year [kgCO2/year]
    CostSewerPerYear [10^4 CNY/year]
    """

    motorPowerInput = 0.00389285714286 * Q * heightDifference  # [kWh/day]  faster
    EnergyUsed = motorPowerInput * 365  # [kWh/year] faster

    PumpGhgPerYear_elec = EnergyUsed * coefficients["GridEmiFactor"]

    # Q_per_s [m3/s]
    Q_per_s = Q/24/3600
    if Q_per_s<= 1:
        PumpCostPerYear = (EnergyUsed * coefficients["Elec_cost"] + coefficients["Pump1"])/10000
        PumpGhgPerYear = PumpGhgPerYear_elec + coefficients["Pump1_ghg"]
    elif Q_per_s<= 3:
        PumpCostPerYear = (EnergyUsed * coefficients["Elec_cost"] + coefficients["Pump2"])/10000
        PumpGhgPerYear = PumpGhgPerYear_elec + coefficients["Pump2_ghg"]
    else:
        PumpCostPerYear = (EnergyUsed * coefficients["Elec_cost"] + coefficients["Pump3"])/10000
        PumpGhgPerYear = PumpGhgPerYear_elec + coefficients["Pump3_ghg"]

    # Error message
    if heightDifference < 0 or PumpGhgPerYear < 0:
        raise Exception("ERROR: Pumping costs cannot be calculated correctly. " + str(heightDifference) + "" + str(
            Q))  # Does not make sense if pumped down

    return PumpGhgPerYear, PumpCostPerYear,



# ===remove
# ===shijingshan_3 and 4 shows that there is no difference between the two getPipeDiameter(Q, slope, stricklerC) methods
def getPipeDiameter(Q, slope, stricklerC):
    """
    This function calculates the pipe diameter according to Manning-Strickler.

    Input Arguments: 
    Q                     -    Flow in pipe
    slope                 -    Slope
    stricklerC            -    Strickler coefficient

    Output Arguments:
    pipeDiameter          -   Needed pipe diameter
    """
    Q = Q/86400.0                                                                           # Convert the flow [m3/day] to [m3/s], 24.0*60.0*60.0  = 86400.0  
    # ===cancel Qmax
    # Qmax = 0.8                                                                              # Maximum filling condition
    normDiameterList = (.25, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.5, 2, 2.5, 3, 4, 6, 8)   # [m] Norm pipe diameters
    # ===add depth_ratio
    depth_ratio = {0.25: 0.55, 0.3: 0.55,
                   0.4: 0.65,
                   0.5: 0.7, 0.6: 0.7, 0.7: 0.7, 0.8: 0.7, 0.9: 0.7,
                   1: 0.75, 1.2: 0.75, 1.5: 0.75, 2: 0.75,
                   2.5: 0.75, 3: 0.75, 4: 0.75, 6: 0.75, 8: 0.75}

    #Iterate list with norm diameters until the calculate flow is bigger
    for D in normDiameterList:

        if slope == 0:  # Error if slope is zero
            # If slope is 0 WWTP is pumped and a Diamter of 0.25 is assumed.
            pipeDiameter = normDiameterList[0]
            return pipeDiameter
            #raise Exception("ERROR: Pipe diamater cannot get calculated because slope is zero.") 
        else:
            # Calculate flow in pipe with diameter D [m3/s]
            #Qfull = stricklerC * (D/4.0)**(2.0/3.0) *math.sqrt(abs(slope)) * (math.pi/4.0) * D**2.0                     # slower
            # cancel Qfull
            # Qfull = stricklerC * (D/4.0)**(0.6666666666666666) *math.sqrt(abs(slope)) * (0.7853981633974483) * D**2.0    # faster

            # central angle of depth
            sitar = math.acos((depth_ratio[D]-0.5)*2)
            # wetted perimeter
            X = (6.28318530718-sitar)*D/2
            # cross-section
            A = (6.28318530718-sitar+math.sin(sitar))*(D**2) / 8
            # hydraulic radius
            R = A/X
            Qfull = A * stricklerC * R**0.6666666666666666 * (abs(slope)**0.5)
            
            # If pipe can bear more than flow as input, select this diameter # 80 % condition
            # ===cancel if Qfull * Qmax >= Q:, change to if Qfull >= Q:
            # if Qfull * Qmax >= Q:
            if Qfull >= Q:
                pipeDiameter = D
                return pipeDiameter
        
        if D == 8: # Not big enough norm-pipe diameter existing
            pipeDiameter = 10
            return pipeDiameter



def getPipeDiameter_Orig(Q, slope, stricklerC):
    """
    This function calculates the pipe diameter according to Manning-Strickler.

    Input Arguments:
    Q                     -    Flow in pipe
    slope                 -    Slope %
    stricklerC            -    Strickler coefficient

    Output Arguments:
    pipeDiameter          -   Needed pipe diameter
    """
    Q = Q / 86400.0  # Convert the flow [m3/day] to [m3/s], 24.0*60.0*60.0  = 86400.0
    # Qmax = 0.8                                                                              # Maximum filling condition
    normDiameterList = (.25, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.5, 2, 2.5, 3, 4, 6, 8)  # [m] Norm pipe diameters
    # ===add depth_ratio
    depth_ratio = {0.25: 0.55, 0.3: 0.55,
                   0.4: 0.65,
                   0.5: 0.7, 0.6: 0.7, 0.7: 0.7, 0.8: 0.7, 0.9: 0.7,
                   1: 0.75, 1.2: 0.75, 1.5: 0.75, 2: 0.75,
                   2.5: 0.75, 3: 0.75, 4: 0.75, 6: 0.75, 8: 0.75}

    # Iterate list with norm diameters until the calculate flow is bigger
    for D in normDiameterList:

        if slope == 0:  # Error if slope is zero
            # If slope is 0 WWTP is pumped and a Diamter of 0.25 is assumed.
            pipeDiameter = normDiameterList[0]
            return pipeDiameter
            # raise Exception("ERROR: Pipe diamater cannot get calculated because slope is zero.")
        else:
            # print(slope)
            # Calculate flow in pipe with diameter D [m3/s]
            Qfull = stricklerC * (D / 4.0) ** (0.6666666666666666) * math.sqrt(abs(slope)) * (
                0.7853981633974483) * D ** 2.0  # faster
            Qmax = depth_ratio[D]

            # If pipe can bear more than flow as input, select this diameter # 80 % condition
            # ===cancel if Qfull * Qmax >= Q:, change to if Qfull >= Q:
            if Qfull * Qmax >= Q:
                pipeDiameter = D
                return pipeDiameter

        if D == 8:  # Not big enough norm-pipe diameter existing
            pipeDiameter = 10
            return pipeDiameter


def costWWTP(flow, coefficients, fc_coefficient_a, fc_coefficient_b):
    """
    This function calculates the costs of a wwtp.
    :param flow: Amount of waste water to be treated [in m3/day]
    :param coefficient_a:  oefficient of scale effect
    :param coefficient_b: oefficient of scale effect
    :param fc_coefficient_a: sensitive analyze of coefficient a
    :param fc_coefficient_b: sensitive analyze of coefficient b
    :return:
    GhgWwtpPerYear - GHG emissions from WWTP per year [kgCO2/year]
    CostWwtpPerYear [10^4CNY/year]
    """
    coefficient_a = coefficients["coefficient_a"]
    coefficient_b = coefficients["coefficient_b"]
    sensFactor_a = (1 + fc_coefficient_a) * coefficient_a
    sensFactor_b = (1 + fc_coefficient_b) * coefficient_b

    # ----GHG----
    # 10^4 ton wastewater per mon
    flow_per_mon = flow * 30 / 10000
    # kg CO2 per t wastewater
    # GhgIntensity = sensFactor_a * (flow_per_mon ** sensFactor_b)
    # kg C02 from WWTP
    # GhgWwtpPerYear = GhgIntensity * flow_per_mon * 12 * 10000
    # combined
    GhgWwtpPerYear = sensFactor_a * (flow_per_mon ** (1 + sensFactor_b)) * 12 * 10000

    # ----CNY----
    # CNY per t wastewater
    # CostIntensity = 2.7538 * (flow ** -0.114)
    # 10^4 CNY from WWTP per year
    # CostWwtpPerYear = CostIntensity * flow_per_mon * 12
    # combined
    CostWwtpPerYear = 5.9635 * (flow ** -0.17) * flow * 365 / 10000

    return GhgWwtpPerYear, CostWwtpPerYear,


def calculateConnectionCosts(totSewerCostI, totPumpCostI,
                             totSewerCostII, totPumpCostII,
                             totSewerCostIII, totPumpCostIII,
                             WWTPcostsI):
    """
    This function calculates the total replacement costs of the system in case a node is connected to the existing system. 
    This is needed in order to compare the costs of a central connection with the reasonable costs.
    
    Input Arguments: 
    expOrMerge                                --    0: In Expansion module, 1: In Expansion module
    pipeCostI, pumpCostIWholePeriodI          --    Pipe-, pumping costs of option I
    pipeCostII, pumpCostIWholePeriodII        --    Pipe- & pumping costs of option II
    pipeCostIII, pumpCostIWholePeriodVIII     --    Pipe- & pumping costs of option III
    
    WWTPcostsI                                --    WWTP costs option I
    
    Output Arguments:
    costConnection                            --    Costs of lowest central connection
    """

    # Calculate total replacement value only of network (including pumps) of the different options
    centralConnectionI = totSewerCostI + totPumpCostI  # total replacement value of pumps and sewer with central connection
    connectionII = totSewerCostII + totPumpCostII  # total replacement value of pumps and sewer without connection
    centralConnectionIII = totSewerCostIII + totPumpCostIII  # total replacement value of pumps and sewer with central connection

    # Calculate connection costs of sewers and pumps not including WWTP. The already existing network needs to be substracted (option II).
    costConnectionI = centralConnectionI - connectionII  # Cost central - cost decentral
    costConnectionIII = centralConnectionIII - connectionII  # Cost central - cost decentral

    # Select lowest connection costs
    if costConnectionI < costConnectionIII:
        costConnection = costConnectionI
    else:
        costConnection = costConnectionIII
    return costConnection


def calculatetotalAnnuities(listWTPs, pumps, sewers, edgeList,
                            nodes, stricklerC, coefficients,
                            fc_SewerCost, fc_coefficient_a, fc_coefficient_b
                            ):
    """
    This function calculates the total system costs of a system
    :param listWTPs: List with wwtps
    :param pumps: List with Pumps
    :param GridEmiFactor: grid emission factor
    :param sewers: Sewers
    :param edgeList: List with edges
    :param nodes: Nodes
    :param stricklerC: Strickler Coefficient
    :param fc_SewerCost: Operation costs
    :param fc_coefficient_a:
    :param fc_coefficient_b:
    :return:
    totSystemCosts      -     Total System Costs
    """

    # calculate WWTPs costs
    completeWWTPGhg = 0
    completeWWTPCosts = 0
    for i in listWTPs:
        GhgWwtp, CostWwtp = costWWTP(
            flow=i[1],
            coefficients=coefficients,
            fc_coefficient_a=fc_coefficient_a,
            fc_coefficient_b=fc_coefficient_b
        )
        completeWWTPGhg += GhgWwtp
        completeWWTPCosts += CostWwtp

    # Calculate pump costs
    completePumpGhg = 0
    completePumpCosts = 0
    for pmp in pumps:
        flow, heightDifference = pmp[1], pmp[2]
        GhgPump, CostPump = getPumpCostsDependingOnFlow(
            Q=flow,
            heightDifference=heightDifference,
            coefficients=coefficients
        )  # pump is found on path
        completePumpGhg += GhgPump
        completePumpCosts += CostPump

    # Calculate sewer costs
    completePipeGhg = 0
    completePipeCosts = 0
    for pipe in sewers:
        if sewers[pipe][0] != ():
            oldNode = pipe
            nextNode = sewers[pipe][0]

            # Get flow
            # for a in nodes:
            #     if a[0] == oldNode:
            #         Q = a[4] + a[8]
            #         break
            a = nodes[oldNode]
            Q = a[4] + a[8]

            # Get distance, slope
            for edge in edgeList:
                if edge[0][0] == oldNode and edge[1][0] == nextNode:  # Stored inverse, thus slope needs to get inverted
                    distance, slope = edge[2], edge[3] * -1  # distance, # slope needs to be inverted
                    break

                if edge[1][0] == oldNode and edge[0][0] == nextNode:
                    distance, slope = edge[2], edge[3]  # distance, # slope stays the same
                    break

            # Get Trench Depth
            # for punkt in nodes:
            #     if punkt[0] == oldNode:
            #         trenchDepthFrom = punkt[3] - punkt[10]
            #         break
            punkt = nodes[oldNode]
            trenchDepthFrom = punkt[3] - punkt[10]

            # for punkt in nodes:
            #     if punkt[0] == nextNode:
            #         trenchDepthTo = punkt[3] - punkt[10]
            #         break
            punkt = nodes[nextNode]
            trenchDepthTo = punkt[3] - punkt[10]


            averageTrenchDepth = (abs(trenchDepthFrom) + abs(trenchDepthTo)) / 2
            # pipeDiameter = getPipeDiameter(Q, slope, stricklerC)
            pipeDiameter = getPipeDiameter_Orig(Q, slope, stricklerC)
            GhgSewer, CostSewer = calculatePipeCosts(
                distance=distance,
                Q=Q,
                pipeDiameter=pipeDiameter,
                slope=slope,
                fc_SewerCost=fc_SewerCost,
                coefficients=coefficients
            )
            completePipeGhg += GhgSewer
            completePipeCosts += CostSewer
    return completePumpGhg, completeWWTPGhg, completePipeGhg, completePumpCosts, completeWWTPCosts, completePipeCosts


def costsPrivateSewers(buildings, buildPoints, pipeDiameterPrivateSewer, averageSlopePrivateSewer, fc_SewerCost, coefficients):
    """
    This function calculates the costs of the private sewers. The private sewers are the closest distance to the street network,
    If the street is too far, the whole distance to the building is used.
    :param buildings: Buildings
    :param buildPoints: Coordinates of Buildings
    :param pipeDiameterPrivateSewer: Pipe Diameter
    :param averageSlopePrivateSewer: Average slope
    :param fc_SewerCost: cost factor sewers
    :return:
    totCostPrivateSewer                -    Total replacement value of private sewers
    """
    costsP_Sewer = 0
    for node in buildings:
        pt_to1_X, pt_to1_Y, gebListe = node[0], node[1], node[2]

        for house in gebListe:
            for geb in buildPoints:
                if geb[0] == house:
                    _, pt_from1_X, pt_from1_Y, Q = geb[0], geb[1], geb[2], geb[4]
                    break

            p0, p1 = (pt_from1_X, pt_from1_Y), (pt_to1_X, pt_to1_Y)
            distance = math.hypot(p0[0] - p1[0], p0[1] - p1[1])

            privateSewercostsPerYear, _ = calculatePipeCosts(
                distance=distance,
                Q=Q,
                pipeDiameter=pipeDiameterPrivateSewer,
                slope=averageSlopePrivateSewer,
                fc_SewerCost=fc_SewerCost,
                coefficients=coefficients,
            )
            costsP_Sewer += privateSewercostsPerYear

    totCostPrivateSewer = costsP_Sewer
    return totCostPrivateSewer


def getCostsOfCrossedWWTPs(allNodesToAddToPN, pathBetweenWWTPs, WWTPS_noCon, sewers_NoCon,
                           nodes_noCon, coefficients,
                           fc_coefficient_a, fc_coefficient_b):
    '''
    This function estimates the costs of all crossed wwtps on the path between two wwtps.
    
    Input:
    allNodesToAddToPN    -    All nodes on the path
    pathBetweenWWTPs     -    Path between WWTP
    WWTPS_noCon          -    WWTP
    sewers_NoCon         -    Sewers
    nodes_noCon          -    Nodes
    EW_Q, coefficient_a, interestRate, fc_wwtpOperation, fc_wwtpReplacement    -    Cost relevant parameters
    
    Output:
    sumCostcrossedWWTP   -    Costs
    '''
    # Iterate path and get the sum of all flow which flows to WWTPs in the path
    allWWTPsInPath = []  # List to store all crossed wwtp with the flow [[ID, flow]]
    sumCostcrossedWWTP = 0  # Total costs

    # Get all WWTPs in Path
    for i in allNodesToAddToPN:
        for wwtp in WWTPS_noCon:
            if i[0] == wwtp[0]:
                if i[0] in pathBetweenWWTPs:
                    allWWTPsInPath.append([i[0], 0])
                    break

    # Iterate path
    if len(allWWTPsInPath) > 0:
        for i in pathBetweenWWTPs:
            wwtpfound = 0

            # get WWTP to which this node flows
            iterate = i
            try:
                while wwtpfound == 0:
                    nextN = sewers_NoCon[iterate]
                    if nextN[0] == ():
                        wwtpfound = 1
                        toWWTP = iterate
                        break
                    iterate = nextN[0]
            except:
                continue  # This node was not in network

            for wwtp in allWWTPsInPath:
                if wwtp[0] == toWWTP:
                    # for n in nodes_noCon:
                    #     if n[0] == i:
                    #         fl = n[8]
                    #         break
                    fl = nodes_noCon[i][8]
                    wwtp[1] += fl
                    break
        for i in allWWTPsInPath:
            costCrossed, _ = costWWTP(
                flow=i[1],
                coefficients=coefficients,
                fc_coefficient_a=fc_coefficient_a,
                fc_coefficient_b=fc_coefficient_b
            )
            sumCostcrossedWWTP += costCrossed
    return sumCostcrossedWWTP
