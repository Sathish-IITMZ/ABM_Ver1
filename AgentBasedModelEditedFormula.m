%% Validation on 15.05.2019 - Predictive Model
% First extracting concentration for every node at every hour.
% For missing hour using simple average interpolation.


%% Reading data - Concentration and Distance
pm25 = readtable(['PM2_5_Hourly.xlsx']);
dist = readtable('Distance22Nodes.xlsx');

pm25(:,1:3) = [];
pm25(1,:) = [];
pm25(:,25) = pm25(:,1);
dist(:,1:3) = [];
dist(1,:) = [];

%% Reading Node Relation Matrix 
load NodeInformationFor22Nodes

%Adjacent Matrix from Node Relation
for i=1:22
   [~,~, relat] = find(NodeRelation(i, 2:5));
    adj(i,relat) = 1;
end
 
[nnodes, ~] = size(adj); %Number of nodes
k = 1;

%Finding the node relation from adjacent matrix;
for i=1:nnodes
   index = find(adj(i,:));
   for j=1: length(index)
       nodeRelation(k,:) = [i index(j)];
       k = k+1;
   end
end

numK = length(nodeRelation); %Number of K's
%% Dij Matrix
DIJ = dist;
for i=1:22
    index = find(adj(i,:)==0);
    DIJ{i,index} = 0;
    index = find(DIJ{i,:});
    DIJ{i,index} = 1./DIJ{i,index}/sum(1./DIJ{i,index});
end
DIJ = table2array(DIJ);

%% Training 

LB = [zeros(numK,1); -100*ones(nnodes,1)];
UB = [ones(numK,1); 100*ones(nnodes,1)];
d = 1;
Pij = zeros(nnodes,24);
Pji = zeros(nnodes,24);
for k=0:23
    bigX = zeros(nnodes,numK);
    bigY = [];
    pm25conc = pm25{:,k+1};
    ConvecMatrix = DIJ.*pm25conc;
    for i=1:nnodes
        index = find(adj(i,:)== 1);
        X = [-ConvecMatrix(i,index), ConvecMatrix(index,i)'];
        [row,col] = find(nodeRelation == i);
        bigX(i,row') = X;
        L = [(pm25{i, d*k+2:d*k+d+1} - pm25{i,d*k+1:d*k+d})]';
        bigY = [bigY; L];
        clear X;
    end   
    bigX = [bigX eye(nnodes)];
    [Par(:,k+1),error(k+1)] = lsqlin(bigX,bigY,[],[],[],[],LB,UB);
    %LB = [zeros(numK,1); (Par(nnodes+1:end,k+1)-10) ];
    %UB = [ones(numK,1); (Par(nnodes+1:end,k+1)+10)];
    for j = 1:nnodes
        [row,col] = find(nodeRelation == j);
        Pij(j,k+1) = bigX(j,row(col==1)')*Par(row(col==1),k+1);
        Pji(j,k+1) = bigX(j,row(col==2)')*Par(row(col==2),k+1);
    end
end
RMSE = sqrt(error/nnodes);

 
% %% Extract the test data - 1st June 2019 
% 
% pmtest25 = zeros(nnodes,24);
% for i=1:22
%     Filename = ['node' num2str(i) '.xlsx'];
%     point = readtable(Filename);
%     for j = 0:23
%        % h = find(point{:,12}.Day >= 22 & point{:,12}.Month >= 5 & ...
%         %  point{:,12}.Hour == j);
%        % Test data on - 22 May
%        h = find(point{:,12}.Day == 1 & point{:,12}.Month == 6 & ...
%            point{:,12}.Hour == j);
%         if ~isempty(h)
%            pmtest25(i,j+1) = mean(point{h,40});
%         end
%     end
%     clear point;
% end
% %pm25(:,25) = pm25(:,1);
% load hourMean;
% %How do we fill in the missing data
% for i=1:22
%     index = find(pmtest25(i,:) == 0);
%     pmtest25(i,index) = hourMean(i,index);
% end
% pm25mean = sum(pmtest25,'all')/(22*24);
% %pmtest25(:,25) = pmtest25(:,1);
% 
% %% Testing-  Model Validation 
% 
% yhat = zeros(nnodes,25);
% yhat(:,1) = pmtest25(:,1);
% 
% for k = 0:22
%     bigX = zeros(nnodes,numK);
%     bigInter = [];
%     for i = 1:22
%        index = find(adj(i,:)== 1);
%         Dij = (1./dist{i,index})/sum(1./dist{i,index});
%         Dlen = length(Dij);  
%         X = -yhat(i,d*k+1:d*k+d)'*Dij;
%            for j=1:Dlen
%                  X = [X  yhat(index(j),d*k+1:d*k+d)'*Dij(j)];
%            end
%            [row,col] = find(nodeRelation ==i);
%            bigX(i,row') = X;
%            Inter = [ yhat(i,d*k+1:d*k+d)]';
%            bigInter = [bigInter;Inter];
%     end
%     bigX = [bigX eye(nnodes)];
%     yhat(:,k+2) = bigX*Par(:,k+1) + bigInter;
% end
% 
% yhat = yhat(:,2:24); %Removing the sensor value, that is initial value.
% pmtest25 = pmtest25(:,2:24);
% 
% for j=1:23
%     RMSE(j) = sqrt((norm(yhat(:,j) - pmtest25(:,j)))^2/22);
%     spear(j) = corr(yhat(:,j),pmtest25(:,j),'Type','Spearman');
%     R2Valid(j) = 1 -(((norm( yhat(:,j)-pmtest25(:,j)))^2)/((norm( pmtest25(:,j)-mean(pmtest25(:,j))))^2));
% end
% %Node metric - 22 nodes
% for j=1:nnodes
%     RMSENode(j) = sqrt((norm(yhat(j,:)' - pmtest25(j,:)'))^2/22);
%     spearNode(j) = corr(yhat(j,:)',pmtest25(j,:)','Type','Spearman');
%     R2ValidNode(j) = 1 -(((norm( yhat(j,:)'-pmtest25(j,:)'))^2)/((norm( pmtest25(j,:)'-mean(pmtest25(j,:)')))^2));
% end
% 
% for i=1:nnodes
%     figure(i)
%     plot(1:23, yhat(i,:),'b','LineWidth',2)
%     hold on
%     plot(1:23,pmtest25(i,:),'--r','LineWidth',2)
%     xlabel('Time')
%     ylabel('PM2.5 Conc')
%     t =['R^2 = ',num2str(R2ValidNode(i))];%,'FontWeight','Bold','FontSize',14);
%     title(i,sprintf(t))
%     grid on
%     legend('Predicted','Observed')
%     ax = gca;
%     set(ax,'FontWeight','bold','LineWidth',2,'FontSize',12,'XLim',[-inf 24])
%     hold off
% end
%%
% ycord = [0;1;2;3;4;4;3;2;4;4;3;2;2;2;1.5;1;0;1;0;0;1;0];
% xcord = [0;0;0;0;0;1;1;1;2;3;2;2;3;4;4;4;4;3;3;2;1;1];
% G = digraph(adj);
% plot(G, 'XData', xcord, 'YData', ycord)
% G.Edges.Weights = Par(nnodes,1);




%%
for jk = 1:24
    ycord = 4*[0;1;2;3;4;4;3;2;4;4;3;2;2;2;1.5;1;0;1;0;0;1;0];
    xcord = 4*[0;0;0;0;0;1;1;1;2;3;2;2;3;4;4;4;4;3;3;2;1;1];
    G = digraph(adj);
    G.Edges.Weight = Par(1:numK,jk); 
for i = 1:22
    for j = 1:22
        ind1 = 0;
        ind2 = 0;
        h = 0;
        k = 0;
        if i ~= j
            ind1 = findedge(G,i,j);
            if ind1 ~= 0
                h = G.Edges.Weight(ind1)*DIJ(i,j)*pm25{i,jk}; 
            end
            ind2 = findedge(G,j,i);
            if ind2 ~= 0
                k = G.Edges.Weight(ind2)*DIJ(j,i)*pm25{j,jk};
            end
            if h ~= 0 && k ~= 0
                if h > k
                    G.Edges.Weight(ind1) = h-k;
                    G = rmedge(G,j,i);
                else 
                    G.Edges.Weight(ind2) = k-h;
                    G = rmedge(G,i,j);                    
                end
            end
        end
    end
end
color = [];
for hk=1:22
    if Par(54+hk,jk)> 0 
        color(hk,:) = [0.8500 0.3250 0.0980];
    else
        color(hk,:) = [0 0.4470 0.7410];
    end
    pmgrad25(hk,jk) = pm25{hk,jk+1} - pm25{hk,jk} ;
    if pm25{hk,jk+1} - pm25{hk,jk} > 0
        color1(hk,:) = [0.8500 0.3250 0.0980];
    else 
        color1(hk,:) = [0 0.4470 0.7410];
    end
end


figure(jk)
h = plot(G,'XData',xcord,'YData',ycord,'EdgeLabel',G.Edges.Weight,'ArrowSize',20,'EdgeFontSize',12,'EdgeFontSize',12,'NodeColor','r','EdgeAlpha',1)
h.NodeFontSize = 12;
title('K Map')
% tiledlayout(1,3,'TileSpacing','none','Padding','tight')
% nexttile
% plot(G,'XData',xcord,'YData',ycord,'MarkerSize',pm25{:,jk},'Marker','s','ShowArrows','off')
% title('PM2.5 Concentration')
% nexttile
% nodeID = [1:22];
% h = plot(G, 'XData', xcord, 'YData', ycord,'MarkerSize',abs(Par(55:end,jk))*3, ...
%      'NodeColor',color,'EdgeFontSize',10,'EdgeColor','k','NodeFontSize',12,'Marker','s','ShowArrows','off')
% labelnode(h,nodeID,cellstr(num2str(Par(55:end,jk))))
% title(BP)
% nexttile
% nodeID = [1:22];
% h = plot(G,'XData', xcord,'YData', ycord,'Marker','s','NodeColor',color1,'EdgeLabel',G.Edges.Weight,'ArrowSize',20,'MarkerSize',abs(pmgrad25(:,jk))*3)
% labelnode(h,nodeID,cellstr(num2str(pmgrad25(:,jk))))
% h.NodeFontSize = 12;
% title('Concentration Difference')
% grid off
end

%% Writing the Estimated Parameter to the excel file

%writematrix(Par,'Parameter.xlsx','Sheet',1)
  


%% Plotting to see result

for i=1:22
    figure(i)
    plot(1:24,pm25{i,2:25},'-+b')
    hold on
    plot(1:24,Par(54+i,1:24),'-*r')
    hold on
    plot(1:24,abs(Pij(i,1:24)),'-x')
    hold on
    plot(1:24,Pji(i,1:24),'-s')
    hold on
    legend('PM2.5','BP','Pij','Pji')
    title(i)
    ylabel('PM2.5 Conc')
    xlabel('Time')
    grid on
    ax = gca;
    set(ax,'FontWeight','bold','LineWidth',2,'FontSize',12,'XLim',[-inf 24])
    hold off
end
%%
% for i=1:nnodes
%     figure(i)
%     plot(yhat(:,i),pmtest25(:,i))
% end
%% Lets do Analysis on BP and Pij,Pji

BP = Par(55:end,:);
BP = BP - mean(BP,2);
Pij_duplicate= Pij - mean(Pij,2);
Pji_duplicate = Pji - mean(Pji,2);
P = pm25{:,2:25};
BPabs = BP./P;
Pij_duplicate_abs = Pij_duplicate./P;
Pji_duplicate_abs = Pji_duplicate./P;
Pt_1 = pm25{:,1:24};
Pt_1 = Pt_1 - mean(Pt_1,2);
Pt_1 = Pt_1./P;
Fraction = [(1:22)' mean(abs(BPabs),2) mean(abs(Pij_duplicate_abs),2) mean(abs(Pji_duplicate_abs),2)];
BPsq = BP.^2./P.^2;
Pij_duplicate_sq = Pij_duplicate.^2./P.^2;
Pji_duplicate_sq = Pji_duplicate.^2./P.^2;
Fractionsq = [(1:22)' mean(BPsq,2) mean(Pij_duplicate_sq,2) mean(Pji_duplicate_sq,2)];
deno =  mean((P - mean(P,2)).^2,2);
Fraction3 = [(1:22)' mean(BP.^2,2)./deno ...
    mean(Pij_duplicate.^2,2)./deno mean(Pji_duplicate.^2,2)./deno];

conv = Pji + Pij;
conv = conv - mean(conv,2);
Fraction4 = [(1:22)' mean(BP.^2,2)./deno mean(conv.^2,2)./deno];
