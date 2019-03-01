function polarizeDriver
% polarizeDriver -- Creates data and plots for analyzing the polarization simulation.
% 
% BY Zoe Koch, 2/28/19

trials = 50;
dim = 4;
n = 20;
maxT = 50;

p_opt = 8;
m_opt = 8;

P = linspace(.501,.8,p_opt);
M = linspace(1,3,m_opt);
% D = (1:2:dim).^2; 
D = 1:dim^2;
d_opt = length(D);

C = zeros(length(D), p_opt, m_opt, trials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN SIMULATIONS

tic

% Display a sample simulation. 
fastpolarize(dim,maxT);

% Get data on the final proportion correct.
pValues = nan;
updateValues = nan;

for d_ind = 1:d_opt
    for p_ind = 1:p_opt
        for m_ind = 1:m_opt
            for trial = 1:trials
                [pValues, updateValues, c] = ...
                    fastpolarize(dim,maxT,D(d_ind),P(p_ind),n,M(m_ind),...
                    'pValues',pValues,'updateValues',updateValues,'display',false);
                C(d_ind, p_ind, m_ind, trial) = c;
            end
%           Reset the saved values when the inputs to fastpolarize change.
            pValues = nan;
            updateValues = nan;
        end
    end
end

% Get data on how the entropy changes through time.
dim = 4;
maxT = 25;
maxD = dim^2;
trials = 1000;
E = zeros(trials,maxT,maxD);
pValues = nan;
updateValues = nan;
for trial = 1:trials
    for d=1:maxD
        [pValues, updateValues, E(trial,:,d)] = ...
            fastpolarize(dim,maxT,d,'display',false,'output','entropy',...
            'pValues',pValues,'updateValues',updateValues);
    end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FIGURES

% FINAL STATE PLOTS

AConsensus = C==0;
BConsensus = C==1;
Polarized = 1 - BConsensus - AConsensus;

% Plot final states vs. D
figure
plot(D,mean(mean(mean(BConsensus,2),3),4),...
    D,mean(mean(mean(AConsensus,2),3),4),...
    D,mean(mean(mean(Polarized,2),3),4))
xlabel('Degree')
ylabel('Average Percentage')
xticks(D)
L = 0:.1:1;
labels = reshape(sprintf('%5.1f%%',L*100),6,[]).';
yticklabels(labels)
suptitle('Final States')
legend('Consensus to the right answer','Consensus to the wrong answer',...
    'Polarization') 
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/statevsD.svg')

% Plot proportion polarized vs. D, P
figure
subplot(1,2,1)
surf(D,P,mean(mean(Polarized,3),4)')
zlabel('Proportion Polarized')
xlabel('Degree')
ylabel({'$P_b$'},'Interpreter','latex')

subplot(1,2,2)
contourf(D,P,mean(mean(Polarized,3),4)')
xlabel('Degree')
ylabel({'$P_b$'},'Interpreter','latex')
colorbar
suptitle('Proportion Ending in a Polarized State')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/polarvsDP.svg')

% Plot proportion polarized vs. D, M
figure
subplot(1,2,1)
surf(D,M,squeeze(mean(mean(Polarized,2),4))')
zlabel('Proportion Polarized')
xlabel('Degree')
ylabel({'$m$'},'Interpreter','latex')

subplot(1,2,2)
contourf(D,M,squeeze(mean(mean(Polarized,2),4))')
xlabel('Degree')
ylabel({'$m$'},'Interpreter','latex')
colorbar
suptitle('Proportion Ending in a Polarized State')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/polarvsDM.svg')


% FINAL PROPORTION CORRECT PLOTS
C_avg = mean(C, 4);

% Plot C vs degree, C vs p, C vs m
figure

subplot(3,1,1)
boxplot(reshape(permute(C,[2 3 4 1]) ,[], d_opt));
xlabel('Degree')
ylabel('Correct')
xticklabels(D)

subplot(3,1,2)
boxplot(reshape(permute(C,[1 3 4 2]) ,[], p_opt));
xlabel({'$P_B$'},'Interpreter','latex')
ylabel('Correct')
xticklabels(round(P,3))

subplot(3,1,3)
boxplot(reshape(permute(C,[1 2 4 3]) ,[], m_opt));
xlabel({'$m$'},'Interpreter','latex')
ylabel('Correct')
xticklabels(round(M,2))

suptitle('Average Final Proportion Correct')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/CvsDvsPvsM.svg')

% Plot C vs. p, m. Degree fixed at its max.
figure
suptitle('Average Final Proportion Correct')

subplot(1,2,1)
surf(P,M,reshape(C_avg(d_opt,:,:),m_opt,p_opt))
xlabel({'$P_B$'},'Interpreter','latex')
ylabel({'$m$'},'Interpreter','latex')

subplot(1,2,2)
contourf(P,M,reshape(C_avg(d_opt,:,:),m_opt,p_opt))
xlabel({'$P_B$'},'Interpreter','latex')
ylabel({'$m$'},'Interpreter','latex')
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/CvsPM.svg')

% Plot C vs. degree, m. P is fixed at ~.8.
figure
suptitle('Average Final Proportion Correct')

subplot(1,2,1)
surf(D,M,reshape(C_avg(:,p_opt,:),m_opt,d_opt))
zlabel('Proportion Correct')
xlabel('Degree')
ylabel({'$m$'},'Interpreter','latex')
xticklabels(D)

subplot(1,2,2)
contourf(D,M,reshape(C_avg(:,p_opt,:),m_opt,d_opt))
zlabel('Proportion Correct')
xlabel('Degree')
ylabel({'$m$'},'Interpreter','latex')
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/CvsDM.svg')

% Plot C vs. degree, p. M is fixed at ~1.
figure
suptitle('Average Final Proportion Correct')

subplot(1,2,1)
surf(1:d_opt,P,reshape(C_avg(:,:,1),p_opt,d_opt))
zlabel('Proportion Correct')
xlabel('Degree')
ylabel({'$P_B$'},'Interpreter','latex')

subplot(1,2,2)
contourf(1:d_opt,P,reshape(C_avg(:,:,1),p_opt,d_opt))
zlabel('Proportion Correct')
xlabel('Degree')
ylabel({'$P_B$'},'Interpreter','latex')
colorbar
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/CvsDP.svg')

% ENTROPY PLOTSa

figure
plot(squeeze(mean(E,1)),'LineWidth',3)
title('Entropy over Time')
xlabel('Time')
ylabel('Entropy')
legend(string(1:maxD),'Location','northeast')
lgd = legend;
lgd.Title.String = 'Degree';
set(findall(gcf,'-property','FontSize'),'FontSize',18)
saveas(gcf,'images/entropyvsT.svg')

