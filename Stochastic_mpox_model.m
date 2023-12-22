%% Calibration of the past mpox outbreak in the Netherlands
% https://www.medrxiv.org/content/10.1101/2023.08.21.23293147v2
% The risk of future mpox outbreaks among men who have sex with men: 
% a modelling study based on cross-sectional seroprevalence data
% Marc C. Shamier, Luca M. Zaeck, Hannelore M. Götz, Bruno Vieyra, 
% Babs E. Verstrepen, Koen Wijnans, Matthijs R.A. Welkers, Elske Hoornenborg, 
% Martin E. van Royen, Kai J. Jonas, Marion P.G. Koopmans, Rory D. de Vries, 
% David A.M.C. van de Vijver, Corine H. GeurtsvanKessel

% the method uses the Gillespie algorithm
% code by David van de Vijver, Viroscience, Erasmus MC, Rotterdam, NL
% d.vandevijver@erasmusmc.nl

clear % clear data in workspace
clc % clear screen
tic % start clock

%% initialize simulations
% set the number of runs
runs = 10^5;

% define matrix used for keeping track of the events
SumMatrix = zeros(10000, 23);

% define event matrix
load EventMatrixModel6.mat
Event = table2array(EventMatrixModel6);

MaxTime = 230; % end time

File = 2;
possible_hit=0; 
%% start iterations
for sim = 1:runs

% set seed used to generate random numbers
rng(sim)

%% initialize data
N= 45000 + rand()*15000; % population at risk
TimeToScaleUpDiagnosis = 7 + rand() *42; % time when PCR testing was available widely
StartWeekOutbreak = 21 - round(TimeToScaleUpDiagnosis/7); % week in which outbreak started
t = 0; % t = time
propR = 0.1 + rand()*0.2; % initial # individuals w/ historical vaccin (13.7% at GGD born before 1974).
seed = 1 + rand()*9; % seeding of an outbreak
E = 0; % # MSM exposed at start outbreak
I = floor(seed); % infected with mpox at start outbreak
R = 0; % immune against mpox due to infection
Vhist = floor(propR * N); % historically vaccinated
Vnew = 0; % newly vaccinated
cumulI = 0; % cumulative number of infections
S = N - E - I - R - Vhist - Vnew; % number of MSM susceptible to mpox
upperBeta = 0 + rand()*0.95; % upperrange beta
upperGamma1= floor(7 + rand()*14); % upperrange time to isolation, period 1
upperGamma2= floor(1 + rand()*7); % upperrange time to isolation, period 2
reduction = 1;

%% data on vaccination
% effectiveness of new vaccine
% https://www.cdc.gov/mmwr/volumes/71/wr/mm7149a5.htm?s_cid=mm7149a5_w
% During this study period, mpox incidence (cases per 100,000 population at risk) 
% among unvaccinated persons was 7.4 (95% CI = 6.0–9.1) times that among 
% persons who received only 1 dose of JYNNEOS vaccine ≥14 days earlier and 
% 9.6 (95% CI = 6.9–13.2) times that among persons who received 
% dose 2 ≥14 days earlier.
% dose 1 1/ 7.4 = 0.1351 (95%CI 0.1099-0.1667)
% dose 2 1/9.6 = 0.1042 (95% CI 0.0758-0.1449) 
VEhist_write = 0.75 + rand()*0.2; % effectiveness historical smallpox vaccination
VEhist = 1 - VEhist_write;
newvaccine = 240 + 120*rand(); % mean daily number vaccinated
% # at least 14 days after vaccination is 78% (95% CI 54 to 89), and during 0-13
% # days after vaccination -4% (-50 to 29) 
VEnew = makedist('Triangular','A',0.475,'B',0.885,'C',0.91); % vaccine efficacy historical smallpox
VEnew1 = 1 - random(VEnew,1,1);
VEnew2 = VEnew1;
VEnew = VEnew1; % no accurate data on second vaccination, therefore take lowest estimate

I7=0; % infected after seven days
I30=0; % infected after 30 days
I60=0; % infected after 60 days
I90=0; % infected after 90 days
I150=0; % infected after 150 days
Ioneweek =0;
Ionemonth =0;
Itwomonth =0;
Iend=0;

%% define initial values matrix that keeps track of epidemic
row = sim - (File - 1) * 10000; 

SumMatrix(row,1) = t; % Column using time
SumMatrix(row,2) = S; % column denoting susceptibles 
SumMatrix(row,3) = Vhist;
SumMatrix(row,4) = Vnew;
SumMatrix(row,5) = E;
SumMatrix(row,6) = I;
SumMatrix(row,7) = R;
SumMatrix(row,8)= cumulI;
SumMatrix(row,9)= I7;
SumMatrix(row,10)= I30;
SumMatrix(row,11)= I60;
SumMatrix(row,12)= I90;
SumMatrix(row,13)= I150;
SumMatrix(row,14)= sim;
SumMatrix(row,15)= N;
SumMatrix(row,16)= propR;
SumMatrix(row,17)= seed;
SumMatrix(row,18)= VEhist;
SumMatrix(row,19)= VEnew1;
SumMatrix(row,20)= upperBeta;
SumMatrix(row,21)= upperGamma1;
SumMatrix(row,22)= upperGamma2;
SumMatrix(row,23)= newvaccine;
SumMatrix(row,24)= StartWeekOutbreak;


ColUpdate = 2:7; % columns including events
ColNewInfect = 6; % column including infection

rates = zeros(1,5); % initialize matrix with rates of events

%% run the Gillespie algoritm
% the code runs until the first of the following events occurred:
% 1) reaching 300 days after start of the epidemic
% 2) all suceptibles are infected
% 3) there are no infected individuals 
% 4) unrealistically high number of infected after 7, 30, 90 or 150 days
while t < MaxTime && S >0 && I > 0 && I7<50 && I30<150 && I90<1550 && I150<1900  
    beta = 0.05 + upperBeta*rand();
    % # serial time
    % # based on Ward et al. BMJ 2022 who reported a serial interval time of 
    % # 8.0 days (95% CI 6.5-9.9) gamma distribution (shape=0.8, rate=0.1)
    alpha = 1/gamrnd(.8,10);
    if t > TimeToScaleUpDiagnosis
        gamma = 1 / (1 + rand() * upperGamma2);
        reduction = 0.5 + 0.5 * rand(); 
        beta = beta * reduction;
    else
        gamma = 1 / (14 + rand() * (upperGamma1-14));
    
    end
    % new infection in unvaccinated individuals
    rates(1) = beta * S * I / N; % by I, S - 1, E + 1
    % new infection in historically vaccinated individual 
    rates(2) = beta * VEhist * Vhist * I  / N; % by I, Vhist - 1, E + 1
    % new infection in newly vaccinated individual 
    rates(3) = beta * VEnew * Vnew * I / N; % by I, Vnew - 1, E + 1
    % disease progression
    rates(4) = alpha * E; % E-1,I+1
    rates(5) = gamma * I; % I-1,R+1

    rate_sum = sum(rates); % Gillespie step 1 calculate sum of rates

    % determine which event is taking place (random event)
    prob_event = zeros(5,4);
    % normalize events to a range {0,1}
    prob_event(:,1) = (rates / rate_sum)';
    % add column with cumulative sum of probabilities
    prob_event(:,2) = cumsum(prob_event(:,1));
    % pick random number to determine which event occurs
    event = rand();
    % find event that occurs
    prob_event(:,3) = prob_event(:,2)>event;
    pickevent = find(prob_event(:,3)>0);
    prob_event(min(pickevent),4) = 1;

    % next line only when statistical toolbox is not available
    % r1 = rand(1); 
    % tau = (1/rate_sum) * log(1/r1);  %time of next event
 
    tau = exprnd(1/rate_sum); % Gillespie step 2 determine when event takes place
    t = t + tau;

    % calculate event that will occur
    WhichEvent = prob_event(:,4) .* Event;
    % Gillespie step 3 determine which event will take place
    UpdateNumbers = sum(WhichEvent,1);
    SumMatrix(row,1) = t;
    % Gillespie step 3 update matrix 
    SumMatrix(row,ColUpdate) = SumMatrix(row,ColUpdate) + ...
        UpdateNumbers(1,ColUpdate);

    % determine if a new infection occurred
    if sum(UpdateNumbers(1,ColNewInfect)) == 1 
        SumMatrix(row,8) = SumMatrix(row,8) + 1;
    end

    % determine number of new infections that have occured at different
    % time points. This is used to speed up code and to avoid that code
    % continues in runs that have an unrealistically high number of
    % infections.
    if t>7 && I7 == 0
        SumMatrix(row,9) = SumMatrix(row,8);
        I7=1;
    end
    if t>30 && I30 == 0
        SumMatrix(row,10) = SumMatrix(row,8);
        I30=1;
    end
    if t>60 && I60 == 0
        SumMatrix(row,11) = SumMatrix(row,8);
        I60=1;
    end
    if t>90 && I90 == 0
        SumMatrix(row,12) = SumMatrix(row,8);
        I90=1;
    end
    if t>150 && I150 == 0
        SumMatrix(row,13) = SumMatrix(row,8);
        I150=1;
    end

    % update new vaccinations
    % outbreak started in week 14 (April 4th), 15 (Apr 11) or 16 (Apr 18)
    % start July 22 (week 29) end Sept 20th (week 38)
    % start period is max 29-14 = 15w (105 d), min 29-16 = 13w (91 d)
    % end period is max 38-14=24w (168d), min 38-16 = 22w (154 d)
    MinVaccDate = (29 - StartWeekOutbreak) * 7;
    MaxVaccDate = (38 - StartWeekOutbreak) * 7;
    if t> MinVaccDate && t<MaxVaccDate
        % average about 125 vaccinations a day 
        SumMatrix(row,4) = SumMatrix(row,4) + tau*newvaccine;
        VaccinS = tau*newvaccine * (1-propR);
        VaccinHist = tau*newvaccine * propR;
        SumMatrix(row,2) = SumMatrix(row,2) - VaccinS;
        SumMatrix(row,3) = SumMatrix(row,3) - VaccinHist;
    end

    % update matrix including all individuals and their mpox/vaccination status 
    SumMatrix(row,1) = t; % Column using time
    S = SumMatrix(row,2);
    Vhist = SumMatrix(row,3);
    Vnew = SumMatrix(row,4);
    E = SumMatrix(row,5);
    I = SumMatrix(row,6);
    R = SumMatrix(row,7);
    Ioneweek = SumMatrix(row,9);
    Ionemonth = SumMatrix(row,10);
    Itwomonth = SumMatrix(row,11);
    beta;

end

% plot(t,cumulI) % uncomment to plot

% small code to put output on screen to see that code is running and does
% not get stuck
if sim == 100
    toc
    SumMatrix;
end

% write output to file
if floor(sim/10000) == sim/10000
    disp(sim)
    FileName = sprintf('Calibration_%d.csv', File);
        File = File +1;
        writematrix(SumMatrix, FileName,'Delimiter',';')
        SumMatrix = zeros(10000, 23);
    toc
end

end
