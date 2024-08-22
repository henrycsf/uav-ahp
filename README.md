# uav-ahp

UAV deployment for 5G emergency networks with preliminary network planning and decision-making via AHP. Developed for the paper named "Dynamic Deployment of UAVs for 5G Emergency Networks Using Multi-Criteria Decision-Making", 
currently under peer review. It also includes a metaheuristic version, for means of comparison. More info in the paper.

The functions drone_positioning_AHP and CS_positioning are the UAV deployment, decision-making part based in AHP and CS, respectively. The file numerology.m defines the network planning and 5G specifications of the network.

To run the codes, select one of the "main" files. "main.m" runs a one-time result of the simulations, selecting different priority cases. "main_plot.m" is made to run simulations several times, take an average of the output metrics,
and put them in bar graphs. "main_comparison.m" is for the complexity and network planning comparison made for Scenario 2 in the paper's results.

Paper Abstract:

UAV Base-Stations (UAV-BS) have proven to be effective in 
emergency networks and high-density events, supported
by 5G and future 6G for enhanced coverage and user demand
management. Traditional deployment optimization of UAVs, an
NP-hard problem, is often too slow for real-time applications,
necessitating faster, adaptable solutions. This paper proposes
transforming the UAV placement problem into a decision-making
process using the Analytic Hierarchy Process (AHP), in an
algorithm named UAV-AHP, where all UAV-BS move along
scanning points in a predetermined path - each with a minimal
distance from one another. Therefore, the proposed UAV-AHP
efficiently suggests optimal UAV-BS positions based on real-time
acquisition of user equipment (UE) signal metrics in scanning
points. Additionally, our study contributes to 5G emergency
network planning, providing solutions for RAN slicing, network
dimensioning and balancing, and automating transmitting power.
Simulations considering three distinct scenarios (rural, urban,
and music festival) are executed in order to validate the proposed
planning and deployment. Results demonstrate that UAV-AHP
outperforms a common bioinspired method for UAV deployment
(Cuckoo Search, or CS) in computation-heavy scenarios, providing 
considerably lower running times and satisfactory solutions
for environments with high density of users.

Index terms: Internet of Flying Things (IoFT), unmanned
aerial vehicle (UAV), emergency network, 5G network planning,
multi-criteria decision-making (MCDM), real-time applications.
