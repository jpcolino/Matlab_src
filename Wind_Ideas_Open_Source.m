function Wind_simultaneous_simulation()

    clc
    clear all
    close all
    warning off all
    format long

    disp(' ')
    disp('-----------------------------------------------------------------------------------')
    disp(' ')
    disp(' Simultaneous Simulation of Wind-Power Production and Spot-Prices ')
    disp(' Author: Jesus Perez Colino')
    disp(' ')
    disp('-----------------------------------------------------------------------------------')
    disp(' ')
    %%% Introduce here your own code to load the calendar, production and prices profiles
    disp(' ')
    disp(' Data loaded sucessfully. Press any key to continue... ')
    pause
    disp('-----------------------------------------------------------------------------------')
    
%   -----------------------------------------------------------------------
%%  1. WIND PRODUCTION: Data Loading and Transformations
%   -----------------------------------------------------------------------
    disp('  ')
    disp(' 1. DAILY WIND PRODUCTION: Transformations, Calibration and Simulation')
    disp('  ')
    
    inputs.dummy    = inputs.data(:,7);
    inputs.year     = nonzeros(inputs.dummy .* inputs.data(:,3));
    inputs.month    = nonzeros(inputs.dummy .* inputs.data(:,2));
    inputs.day      = nonzeros(inputs.dummy .* inputs.data(:,1));

    inputs.mwh      = nonzeros(inputs.dummy .* inputs.data(:,6));
    inputs.log_mwh  =  [0; diff(log(inputs.mwh))];
    
    numdata         = size(inputs.mwh,1);
        
    % Spliting wind production by years
    
    years = min(inputs.year) + 1 : max(inputs.year) - 1;
    numyears = size(years,2);
    
    valueset = zeros( numdata , numyears );
    
    count = 0;
    for year = years(1,1) : years(1,end)
        clear n m i j
        count = count + 1; 
        [n , m] = find ( inputs.year == year );
        for i = n(1) : numdata
            if inputs.year(i,1) == year
               valueset(i,count) = inputs.log_mwh(i,1);
            else
                break
            end
        end
    end
    numyd = 10000000;
    for i = 1:numyears;
        eval(['inputs_log_',num2str(years(1,i)),'= nonzeros(valueset(:,i));']);
        eval(['numyd = min(numyd,size(inputs_log_',num2str(years(1,i)),',1));']);
    end
    for i = 1:numyears;
        eval(['inputs_log_',num2str(years(1,i)),...
            '=inputs_log_',num2str(years(1,i)),'(1:numyd,1);']);
    end

   clear valueset         
   
   year_matrixdata = [zeros(1:numyd,numyears) ] ;
   for i = 1:numyears;
        eval(['year_matrixdata(1:numyd,i)=inputs_log_',num2str(years(1,i)),'(1:numyd,1);']);
   end
%   -----------------------------------------------------------------------
%% 2. WIND PRODUCTION: 
%     Non-parametric regressions using Local Regression Smoothing for the SEASONAL COMPONENT
%   -----------------------------------------------------------------------

    disp('-----------------------------------------------------------------------------------')
    disp('  ')
    disp(' 2. WIND PRODUCTION: Non-parametric regressions using Local Regression Smoothing for the SEASONAL COMPONENT')
    disp('  ')
    
    mean_year = mean(year_matrixdata,2);
    
    n = find(inputs.year == 2012 & inputs.month == 1 &...
        inputs.day==1 );
    
    initial_point = inputs.mwh(n,1); 
    
    prod = [    initial_point.*exp(mean_year(1,1))
                exp(mean_year(2:end,1)) ];
            
    max_point = max(inputs.data(:,end));
    
    mean_prod = cumprod(prod);
    time = [1:size(mean_prod,1)]';  
%     curve1 = smooth(mean_prod(:,1),0.75,'rloess') ;   
    curve1 = smooth(mean_prod(:,1),0.50,'rloess') ;  
    curve2 = smooth(mean_prod(:,1),0.30,'rloess') ;  
    curve3 = smooth(mean_prod(:,1),0.10,'rloess') ;  
%     curve2 = smooth(mean_prod(:,1),24*7*4*6/numyd,'rloess') ;
%     curve3 = smooth(mean_prod(:,1),24*7*4*3/numyd,'rloess') ;
%     curve4 = smooth(mean_prod(:,1),24*7,'moving') ;
   close all
   
   figure(1)
   plot(mean_prod,'--b')
   hold on
   plot(curve1,'-r','LineWidth',2)

   hold on
   plot(curve2,'-g','LineWidth',2)
   hold on
   plot(curve3,'-k','LineWidth',2)
%    hold on
%    plot(curve4,'-c','Linewidth',2); grid on;
   title('Diffenrent Smoothing Functions fitted to hourly production')
%    legend('rloess 0.75', 'rloess 24*7*4/numyd','rloess 24*7/numyd')
%   -----------------------------------------------------------------------   
%%  3. WIND PRODUCTION: Hourly Scheme by Blocks/Months
%   -----------------------------------------------------------------------   
%   Nothing here. It does not apply in the daily scheme. 
%   -----------------------------------------------------------------------
%%  4. WIND PRODUCTION: Construction of the Output Vector: Seasonalize drift
%   -----------------------------------------------------------------------

    disp('-----------------------------------------------------------------------------------')
    disp('  ')
    disp(' 4. WIND PRODUCTION: Construction of the Output Vector: Seasonalize drift ')
    disp('  ')


    output.year     = datevec(datenum(2012,01,01):datenum(2013,01,01),'yyyy');
    output.month    = datevec(datenum(2012,01,01):datenum(2013,01,01),'mm');
    output.day      = datevec(datenum(2012,01,01):datenum(2013,01,01),'dd');
    output.seasonal = curve1 ;
    
    figure(3)
    plot(mean_prod,'--b'); grid on;
    hold on
    plot(curve1,'-r','LineWidth',2)
%     hold on
%     plot(output.seasonal ,'-r','LineWidth',1); grid on
    title('Total Seasonal Component')
    
    disp('Work done. Check figure(3)')
%  ------------------------------------------------------------------------ 
%% 5. Calibration of a Log-Normal Orstein-Uhlenbeck process
%  ------------------------------------------------------------------------ 
    disp('-----------------------------------------------------------------------------------')
    disp('  ')
    disp(' 5. Calibration of a Log-Normal Orstein-Uhlenbeck process for WIND PRODUCTION ')
    disp('  ')

    clear x y b stdev revrate resid dummy rnd m n simul diff_l_input   
    final           = output.seasonal;     % <- Seasonal Funtion
    m_input         = mean_prod;           % <- Average Production Projected (Recomended)
    % input        = inputs.mwh;        % <- Full Historical Production(vol = 1450%) 
%  ------------------------------------------------------------------------ 
%  FIRST calibration method: Using least squares regression.
%  That idea is correct based in the transformation to the log.
    y = log(m_input(2:end,end));
    x = [ones(size(m_input,1)-1,1) log(m_input(1:end-1,end))];
    [b,~,resid]   = regress(y,x);
    revrate_1     = -( log(b(2)) );
    stdev_1       = std(resid) * sqrt( -2*log(b(2)) / ( 1-b(2)*b(2) ) );
%  ------------------------------------------------------------------------ 
%  SECOND calibration method: 
%  Standard deviation estimation: Usual estimation
    size_simul   = size(final,1);
    diff_l_input = diff(log(m_input));
    stdev        = sqrt((1/(size(m_input,1)-1)).* sum((diff_l_input).^2));
    % Reversion Rate estimation: 
    revrate_num = 0;
    for i = 2 : size_simul
        revrate_num = revrate_num + ...
            ((log(final(i-1,1))-log(m_input(end-size_simul+i-1,end))).*...
            (log(m_input(end-size_simul+i,end))-log(final(i,1)))/stdev);
    end
    revrate_den = 0;
    for i = 1 : size_simul-1
        revrate_den = revrate_den + ...
            ((log(final(i,1))-log(m_input(end-size_simul+i,end))).*...
            (log(m_input(end-size_simul+i,end))-log(final(i,1)))/stdev);
    end
    revrate = -log(revrate_num/revrate_den);    
%  ------------------------------------------------------------------------ 
%  THIRD Calibration method: JUST INDICATIVE, NOT CONCLUSIVE due to LOG-N
%  Using C.Ball, W.N.Torous (1983) MLE method for a MRJD (external function)

    [alpha_3,beta_3,sigma_3,mu_3,gamma_3,lambda_3] = mrjd_mle(log(m_input));
    revlevel_3   = alpha_3/beta_3;
    revrate_3    = beta_3;
%  ------------------------------------------------------------------------ 
%  Calibration Control LOOP

    if abs((stdev_1 - stdev)/stdev)*100 > 5 || ...
            abs((revrate_1 - revrate)/revrate)*100 > 5
        
       formatSpec1 = 'Stdev is in Calibration 1: %1.4f Calibration 2: %1.4f Calibration 3: %1.4f \n';
       formatSpec2 = 'RevRate is in Calibration 1: %1.4f Calibration 2: %1.4f Calibration 3: %1.4f \n';
       disp('-------------------------------------------------------------------------------------')
       fprintf(formatSpec1,stdev_1,stdev,sigma_3)
       fprintf(formatSpec2,revrate_1, revrate, revrate_3)
       disp('-------------------------------------------------------------------------------------')
       disp('Fail in the calibration: STOP and CHECK')
       disp('or press any key to continue...')
       disp('-------------------------------------------------------------------------------------')
      % pause
    end
    
%  ------------------------------------------------------------------------
%% 6. Simulation of WIND PRODUCTION as a log-normal Ornstein-Uhlenbeck process
% -------------------------------------------------------------------------
    disp('-----------------------------------------------------------------------------------')
    disp('  ')
    disp(' 6. WIND SIMULATION: Simulating a log-OU process with seasonality')
    disp('  ')
    disp(' Work in progress. Please, be patient...' )


    firstdate   = '01-Jan-2013 00:00:00';
    lastdate    = '30-Jun-2028 00:00:00';
    simul_date  = [datenum(firstdate):datenum(lastdate)]';
    simulp.date = datevec(simul_date,'yyyy');
    m = 1000 ;                  % number of simulations
    n = size(simul_date,1) ;     % number of steps

    % final = [repmat
    % m = 1000 ;               % number of simulations
    % n = size(final,1) ;     % number of steps

    clear x; x = [repmat(final,size(simulp.date,1)/365,1)];
    final=[repmat(final,size(simulp.date,1)/365,1); final(1:(size(simulp.date,1)-size(x,1)),1)];


    simul = zeros(size(final,1),m);
    simul(1,:) = log(final(1,1));
    rnd = randn(n,m);
    grad = gradient(log(final));

for i = 2 : size(final,1)
    % simul(i,:) = (revrate.*(log(final(i,1))-(simul(i-1,:)))+stdev.*rnd(i,:));
    % simul(i,:) = grad(i,:) + revrate.*(log(final(i,:))-simul(i-1,:))+stdev.*rnd(i,:);
    % simul(i,:) = revrate.*(grad(i,:) + log(final(i,:))-simul(i-1,:))+stdev.*rnd(i,:);
    simul(i,:) = grad(i,:)+revrate.*(log(final(i,:))-simul(i-1,:))+stdev.*rnd(i,:);
    simul(i,:) = simul(i,:)+simul(i-1,:);
end

    final_simul = exp(simul);
    final_simul = final_simul ./ exp((stdev.^2)/4);
    maxprod = max(m_input(:,end));

for i = 1 : size(final_simul,1)
    for j = 1 : size(final_simul,2)
        if final_simul(i,j) > maxprod
            final_simul(i,j) = maxprod;
        end
    end
end


% [r,c]=find(final_simul>max(input(:,end)));
% final_simul(r,c)=max(input(:,end));

figure(4)

    subplot(1,2,1)
    plot(final_simul(:,1:5))
    title('Simulations vs. Historical data')
    hold on
    % plot(input(end-size(final_simul,1)+1:end,end),'-g','LineWidth',3)
%     plot(m_input(1:size(final_simul,1),end),'-g','LineWidth',2)
    grid on
    
    subplot(1,2,2)
    plot(mean(final_simul,2),'-b')
    hold on
    plot(final,'-r','LineWidth',3)
    hold on
    plot(final,'--r','LineWidth',1)
    title('Drift and Mean of Wind Production')
    grid on

% figure (5)
%     for i = 1 : size(inputs.data,1)
%         if inputs.year(i,1) == 2011
%            inputs.y2011(i,1) = inputs.data(i,end);
%         end
%     end
%     plot(nonzeros(inputs.y2011),'-b','LineWidth',2)
%     hold on
%     plot(output.seasonal ,'-r','LineWidth',1); grid on
    disp('  ')
    disp(' Work done. Please, check figure 4 ' )
    clear input_*
%  ------------------------------------------------------------------------
%% 7. Calculating sensitivities of spot prices against wind production
%  ------------------------------------------------------------------------
    disp('-----------------------------------------------------------------------------------')
    disp('  ')
    disp(' 7. WIND PRODUCTION and SPOT PRICES: Calculating non-linear sensitivities')
    disp('  ')
    disp(' Work in progress. Please, be patient...' )

    % input1 = xlsread('C:\Users\jepco\Dropbox\Wind Project\2012FundamentalUK.xlsx', 'detalj-avregn1', 'C8:Y146'); 
    % input1 = xlsread('C:\Users\suso\Dropbox\Wind Project\2012FundamentalUK.xlsx', 'detalj-avregn1', 'C8:Y146'); 
    wind = input1(:,14);
    nowind = [input1(:,1:13) input1(:,15:end-2)];
    prices = input1(:,end);
    demand = input1(:,end-1);

    % Price and Power production transformartions: relative variations
    ldprices = diff(log(prices));
    ldnowind = diff(log(sum(nowind,2)));
    ldwind   = diff(log(wind));
    lddemand = diff(log(demand));

    % Price and Power production transformartions: absolute variations
    dprices = diff(prices);
    dnowind = diff(sum(nowind,2));
    dwind   = diff(wind);
    ddemand = diff(demand);
    
% Scattering first results about prices and production
% -------------------------------------------------------------------------

figure(6)

    subplot(2,2,1)
    scatter(prices, demand, 'filled')
    title('Scatter price and demand')
    grid on

    subplot(2,2,2)
    scatter(ldprices, lddemand, 'filled'); grid on; hold on;
    brod = robustfit(lddemand, ldprices);
    plot(lddemand,brod(1)+brod(2)*lddemand,'g','LineWidth',2)
    bls = regress(ldprices, [ones(size(lddemand,1),1) lddemand]); hold on;
    plot(lddemand, bls(1)+bls(2)*lddemand,'r','LineWidth',2);
    title('Scatter relative varations of price and production')
    grid on

    subplot(2,2,3)
    clear idx dm mean Cov robust outliers x y
    %scatter(ldnowind, ldprices , 'filled' ); grid on; hold on;
    [idx,~,~,~] = kur_rce([ldnowind ldprices]);
    outliers = [idx idx].*[ldnowind ldprices];
    [x,y] = find(idx == 0);
    robust = zeros(size(x,1), 2);
    for i = 1 : size(x,1)
        robust(i,:) = [ldnowind(x(i),1) ldprices(x(i),1)];
    end
    scatter(robust(:,1), robust(:,2), 'filled' ); grid on; hold on;
    brod = robustfit(ldnowind, ldprices);
    plot(ldnowind,brod(1)+brod(2)*ldnowind,'g','LineWidth',2); hold on;
    bls = regress(robust(:,2),[ones(size(robust(:,1),1),1) robust(:,1)]); 
    plot(robust(:,1),bls(1)+bls(2)*robust(:,1),'r','LineWidth',2); hold on;
    scatter(outliers(:,1), outliers(:,2),'or');
    title('Scatter relative var. price and No-wind')

    subplot(2,2,4)
    clear idx dm mean Cov robust outliers x y d_YFIT X1FIT X2FIT x y
    [idx,~,~,~] = kur_rce([ldwind ldprices]);
    outliers = [idx idx].*[ldwind ldprices];
    [x,y] = find(idx == 0);
    robust = zeros(size(x,1), 3);
    for i = 1 : size(x,1)
        robust(i,:) = [ldwind(x(i),1) ldprices(x(i),1) ldnowind(x(i),1)];
    end
    scatter(robust(:,1), robust(:,2), 'filled' ); grid on; hold on;
    brod = robustfit(robust(:,1), robust(:,2));
    plot(robust(:,1),brod(1)+brod(2)*robust(:,1),'g','LineWidth',2); hold on;
    bls = regress(robust(:,2),[ones(size(robust(:,1),1),1) robust(:,1)]); 
    plot(robust(:,1),bls(1)+bls(2)*robust(:,1),'r','LineWidth',2); hold on;
    scatter(outliers(:,1), outliers(:,2),'or')
    title('Scatter relative var. price and Wind')
    grid on

figure(7)
    % subplot(1,2,1)
    clear idx dm mean Cov robust outliers x y d_YFIT X1FIT X2FIT x y
    [idx,~,~,~] = kur_rce([ldwind ldnowind ldprices]);
    outliers = [idx idx idx].*[ldwind ldnowind ldprices];
    [x,y] = find(idx == 0);
    robust = zeros(size(x,1), 3);
    for i = 1 : size(x,1)
        robust(i,:) = [ldwind(x(i),1) ldnowind(x(i),1) ldprices(x(i),1) ];
    end
    r_wind   = robust(:,1);
    r_nowind = robust(:,2);
    r_prices = robust(:,3);
    scatter3(r_wind , r_nowind, r_prices, 'ob','filled'); hold on;

    % Regression of degree 3 polinomial over the Relative Variations
    
    % X = [ones(size(r_wind)) r_wind r_nowind r_wind.^2 r_wind.^3 r_wind.*r_nowind (r_wind.^2).*r_nowind];
    % X = [ones(size(ldwind)) ldwind ldnowind ldwind.^2 ldwind.^3 ldwind.*ldnowind (ldwind.^2).*ldnowind];
    X = [ones(size(ldwind)) ldwind ldnowind ldwind.^2 ldwind.^3 ldnowind.^2 ldnowind.^3 ldwind.*ldnowind (ldwind.^2).*ldnowind];
    [b,~,~,~,stats]  = regress(ldprices,X);
    x1fit       = min(ldwind) : 1/100 : max(ldwind);
    x2fit       = max(ldnowind):(-(max(ldnowind)-min(ldnowind))/size(x1fit,2)):min(ldnowind);
    % x2fit       = min(ldnowind) : 1/100 : max(ldnowind);
    % x2fit       = -0.05 : 1/100 : 0.05;
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);

    % YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.^2 + b(5)*X1FIT.^3 + b(6)*X1FIT.*X2FIT + b(7)*(X1FIT.^2).*X2FIT  ;
    YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.^2 + b(5)*X1FIT.^3 +...
           b(6).*X2FIT.^2 + b(7).*X2FIT.^3 + b(8).*X1FIT.*X2FIT + b(9)*(X1FIT.^2).*X2FIT  ;

    str = sprintf('Coefficient of Determination %2.4f ', stats(1,1));

    meshc(X1FIT,X2FIT,YFIT); hold on;
    scatter3(outliers(:,1),outliers(:,2),outliers(:,3),'ob')
    text(1,-0.15,0.15,str )
    xlabel('r Wind Production (rv)')
    ylabel('r No Wind Production (rv) (Vertical load)')
    zlabel('r Power Prices (rv)')
    title ('What is linear in the Price-Production relationship')
    view(50,10)
    

    
    % Regression of degree 3 polinomial over the Absolute Variations
    
    % X = [ones(size(r_wind)) r_wind r_nowind r_wind.^2 r_wind.^3 r_wind.*r_nowind (r_wind.^2).*r_nowind];
    % X = [ones(size(ldwind)) ldwind ldnowind ldwind.^2 ldwind.^3 ldwind.*ldnowind (ldwind.^2).*ldnowind];
    X = [ones(size(dwind)) dwind dnowind dwind.^2 dwind.^3 dnowind.^2 dnowind.^3 dwind.*dnowind (dwind.^2).*dnowind];
    [bb,~,~,~,stats]  = regress(dprices,X);
    x1fit       = min(dwind) : 1: max(dwind);
    x2fit       = max(dnowind):(-(max(dnowind)-min(dnowind))/size(x1fit,2)):min(dnowind);
    % x2fit       = min(ldnowind) : 1/100 : max(ldnowind);
    % x2fit       = -0.05 : 1/100 : 0.05;
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);

    % YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.^2 + b(5)*X1FIT.^3 + b(6)*X1FIT.*X2FIT + b(7)*(X1FIT.^2).*X2FIT  ;
    YFIT = bb(1) + bb(2)*X1FIT + bb(3)*X2FIT + bb(4)*X1FIT.^2 + bb(5)*X1FIT.^3 +...
           bb(6).*X2FIT.^2 + bb(7).*X2FIT.^3 + bb(8).*X1FIT.*X2FIT + bb(9)*(X1FIT.^2).*X2FIT  ;

    str = sprintf('Coefficient of Determination %2.4f ', stats(1,1));

    meshc(X1FIT,X2FIT,YFIT); hold on;
    scatter3(outliers(:,1),outliers(:,2),outliers(:,3),'ob')
    text(1,-0.15,0.15,str )
    xlabel('r Wind Production (rv)')
    ylabel('r No Wind Production (rv) (Vertical load)')
    zlabel('r Power Prices (rv)')
    title ('What is linear in the Price-Production relationship')
    view(50,10)
    

    
% Marginal function estimation

    clear d_YFIT X1FIT X2FIT x y
    % x1fit       = min(ldwind) : 1/100 : max(ldwind);
    % x2fit       = 0.04 : -.08/size(x1fit,2) : -0.04;
    x1fit       = 2*min(ldwind):2*(max(ldwind)-min(ldwind))/100:2*max(ldwind);
    % x2fit       = zeros(size(x1fit,1),size(x1fit,2));
    % x2fit       = 2*min(ldnowind):2*((max(ldnowind)-min(ldnowind))/size(x1fit,2)):2*max(ldnowind);
    x2fit       = 3*max(ldnowind):3*(-(max(ldnowind)-min(ldnowind))/size(x1fit,2)):3*min(ldnowind);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    for i = 1 : size(x1fit,2)
        d_YFIT(i,1) = b(1) + b(2).*x1fit(1,i) + b(4).*x1fit(1,i).^2 + b(5).* x1fit(1,i).^3 ...
            + b(8).*x1fit(1,i).*x2fit(1,i)+ b(9)*(x1fit(1,i).^2).*x2fit(1,i) ...
            + b(3)*x2fit(1,i) + b(6).*x2fit(1,i).^2 + b(7).*x2fit(1,i).^3  ;
    %     d_YFIT(i,1) = b(2) + 2 .* b(4) .* x1fit(1,i) + 3 .* b(5) .* x1fit(1,i).^2 + b(8) .*x2fit(1,i)...
    %        + 2.* b(9).*(x1fit(1,i)).*x2fit(1,i)  ;
    end
    
    
figure(8)
subplot(1,2,1)
    plot(X1FIT(1,1:size(d_YFIT,1)),d_YFIT(:,1),'-r','LineWidth',2); grid on;
    xlabel('r Wind Production (rv)')
    % ylabel('r No Wind Production (rv) (Vertical load)')
    ylabel('r Power Prices (rv)')
    title({'First derivative with respect the wind production';'First derivative of the polynomial regression'})    
    
    
    for i = 1 : size(x1fit,2)
        d_YFIT(i,1) = bb(1) + bb(2).*x1fit(1,i) + bb(4).*x1fit(1,i).^2 + bb(5).* x1fit(1,i).^3 ...
            + bb(8).*x1fit(1,i).*x2fit(1,i)+ bb(9)*(x1fit(1,i).^2).*x2fit(1,i) ...
            + bb(3)*x2fit(1,i) + bb(6).*x2fit(1,i).^2 + bb(7).*x2fit(1,i).^3  ;
    end
    
subplot(1,2,2)
    plot(X1FIT(1,1:size(d_YFIT,1)),d_YFIT(:,1),'-r','LineWidth',2); grid on;
    xlabel('r Wind Production (rv)')
    % ylabel('r No Wind Production (rv) (Vertical load)')
    ylabel('r Power Prices (rv)')
    title({'First derivative with respect the wind production';'First derivative of the polynomial regression'})    
    
   disp('Work done! Please check the figures 6, 7 and 8')
   disp(' ')
%  ------------------------------------------------------------------------
%% 8. Spot PRICE Modelling: Calibration and Join simulation 
%  ------------------------------------------------------------------------
    disp('-----------------------------------------------------------------------------------')
    disp(' ')
    disp(' 8. UK Spot PRICE Modelling: calibration and join simulation ')
    disp('  ')
    disp (' Work in progress... please, be patient')
    
    
%%  8.1. Preparing the time series
%  ------------------------------------------------------------------------
    % input2.data is the Daily Spot Prices for an specific market (Weighted average prices)
    
    % Week - DAYS:
    input2.year             = input2.data(:,1) ;
    input2.month            = input2.data(:,2) ;
    input2.day              = input2.data(:,3) ;
    input2.prices           = input2.data(:,end) ; 
    
    % Week - ENDS:
    input2.year_wend         = input2.data_wend(:,1) ;
    input2.month_wend        = input2.data_wend(:,2) ;
    input2.day_wend          = input2.data_wend(:,3) ;
    input2.prices_wend       = input2.data_wend(:,end);
    

    firstdate_wd        = datenum([num2str(input2.day(end,1)) '-'...
                            num2str(input2.month(end,1)) '-'  num2str(input2.year(end,1))],'dd-mm-yyyy');
    firstdate_wend      = datenum([num2str(input2.day_wend(end,1)) '-'...
                            num2str(input2.month_wend(end,1)) '-'  num2str(input2.year_wend(end,1))],'dd-mm-yyyy');
    
    lastdate_wd         = datenum([num2str(input2.day(1,1)) '-'...
                            num2str(input2.month(1,1)) '-'  num2str(input2.year(1,1))], 'dd-mm-yyyy');
    lastdate_wend       = datenum([num2str(input2.day_wend(1,1)) '-'...
                            num2str(input2.month_wend(1,1)) '-'  num2str(input2.year_wend(1,1))], 'dd-mm-yyyy');
                                        
                        
    firstdate           = max (firstdate_wd, firstdate_wend);
    lastdate            = min (lastdate_wd, lastdate_wend);
    hist_data           = [ firstdate : lastdate ]';
    input2.hist_dat     = datevec(hist_data, 'yyyy');
    size_full           = size(input2.hist_dat, 1);
    input2.full_prices  = zeros(size_full, 1);
    
    
    for i = 1 : size_full
        
        date = [num2str(input2.hist_dat(i,1)),'-',num2str(input2.hist_dat(i,2)),'-', num2str(input2.hist_dat(i,3))];

        if weekday(date) == 1 || weekday(date) == 7
  
               for j = 1:size(input2.data_wend,1)
                   if ( input2.data_wend(j,1) == input2.hist_dat(i,1) && ...
                        input2.data_wend(j,2) == input2.hist_dat(i,2) && ...
                        input2.data_wend(j,3) == input2.hist_dat(i,3) );
                       break               
                   end
                   
               end
           
          input2.hist_dat(i,4) = input2.prices_wend(j,1);

        else
           
               for j = 1:size(input2.data,1)
                   if (input2.data(j,1) == input2.hist_dat(i,1) && ...
                       input2.data(j,2) == input2.hist_dat(i,2) && ...
                       input2.data(j,3) == input2.hist_dat(i,3) );
                       break
                   end
               end
            
          input2.hist_dat(i,4) = input2.prices(j,1);
           if j == size(input2.data,1) 
             input2.hist_dat(i,4) = input2.hist_dat(i-1,4); 
           end
                          
           end
    end
     
    numdata             = size(input2.hist_dat,1);
    
    % Splitting week days from week-end
    
    for i = 1 : numdata
        date = [num2str(input2.hist_dat(i,1)),'-',num2str(input2.hist_dat(i,2)),'-', num2str(input2.hist_dat(i,3))];
        if weekday(date) == 1 || weekday(date) == 7
           input2.weekend(i,1) = 1 ;
           input2.weekday(i,1) = 0 ;
        else
           input2.weekend(i,1) = 0 ;
           input2.weekday(i,1) = 1 ;
        end
    end

   %  input2.p_wend   =  input2.weekend .* [0; diff(input2.hist_dat(:,4))];
    input2.p_wd     =  [0; input2.hist_dat(:,4)];
    

    
%%  8.2. Calibration of a MRJD sde for spot hourly power prices
%       dX = (alpha - beta*X)*dt + sigma*dB + N(mu,gamma)*dN(lambda)     
%  ------------------------------------------------------------------------
    
    % Calibration of the Week-day Model using Ball-Torous MLE method
    
    [alpha_p_wd,beta_p_wd,sigma_p_wd,mu_p_wd,gamma_p_wd,lambda_p_wd] = mrjd_mle(input2.hist_dat(:,4));
    revlevel_p_wd    = alpha_p_wd / beta_p_wd;
    revrate_p_wd     = beta_p_wd;
    
    %  Alternative calibration method for prices: Using least squares regression.
    clear y x f
    y       = input2.hist_dat(2:end,4);
    x       = [ones(size(input2.hist_dat(1:end-1,4),1),1) (input2.hist_dat(1:end-1,4))];
    [f,~,resid]   = regress(y,x);
    revrate_1_p_wd     = -( log(f(2)) );
    revlevel1_p_wd     = f(1)/(1-f(2));
    stdev_1_p_wd       = std(resid) * sqrt( -2*log(f(2)) / ( 1-f(2)*f(2) ) );
    
    
    
    
    % Calibration of Week-end Model
    
%     [alpha_wend,beta_wend,sigma_wend,mu_wend,gamma_wend,lambda_wend] = mrjd_mle_d(nonzeros(input2.p_wend));
%     revlevel_wend    = alpha_wend/beta_wend;
%     revrate_wend     = beta_wend;
    
    
% 8.3. Estimating the Price Impact from the UK wind production     
%  ------------------------------------------------------------------------

    clear d_YFIT X1FIT X2FIT x y
    x1fit       = diff(final_simul);
    x2fit       = zeros(size(x1fit,1),size(x1fit,2));
    %x2fit       = 0;
    
    for j = 1 : size(x1fit,2);
        for i = 1 : size(x1fit,1)
            d_YFIT(i,j) = bb(1) + bb(2).*x1fit(i,j) + bb(4).*x1fit(i,j).^2 + bb(5).* x1fit(i,j).^3 ...
                        + bb(8).*x1fit(i,j).*x2fit(i,j) + bb(9)*(x1fit(i,j).^2).*x2fit(i,j) ...
                        + bb(3)*x2fit(i,j) + bb(6).*x2fit(i,j).^2 + bb(7).*x2fit(i,j).^3 + 0.23;
        end
    end
% 
% figure(9)
%     subplot(3,1,1)
%     plot(final_simul(1:size(d_YFIT,1),1),'-g','LineWidth',2); grid on; 
%     xlabel('Hours in a year')
%     ylabel('Variation of Wind Production (% of Mwh)')
%     title({'Simulation of hourly Wind Production (in green)'})
% 
%     subplot(3,1,2)
%     plot(x1fit(1:size(d_YFIT,1),1),'-b','LineWidth',1); grid on; 
%     xlabel('Hours in a year')
%     ylabel('Variation of Wind Production (% of Mwh)')
%     title({'Simulation of hourly variation of Wind Production (in blue)'})
%     
%     subplot(3,1,3)
%     plot(d_YFIT(:,1),'-r','LineWidth',2); grid on
%     xlabel('Hours in year')
%     ylabel('Jump of Spot Power Prices (GBP/MWH)')
%     title({'Simulation of impact in Spot Prices produced by changes in Wind (in red)'})

% 8.4. Simulating the continous part of the Spot Prices, with three
%      different frameworks: Peak, Off-peak and Weekends
% ------------------------------------------------------------------------
    firstdate   = '01-Jan-2013 ';
    % firstdate   = today();
    lastdate    = '30-Jun-2028 ';
    simul_date  = [datenum(firstdate):datenum(lastdate)]';
    simulp.date = datevec(simul_date,'yyyy');
    m = 1000 ;              % number of simulations
    
   %  n = min(size(simul_date,1),size(final_simul,1)) ;     % number of steps
    n = size(simul_date,1) ; 
    simulp.date = simulp.date(1:n,:) ;
    
    rnd = randn(n,m);
    
    
%     simulp.hour     = simulp.date(:,4);
    simulp.day      = simulp.date(:,3);
    simulp.month    = simulp.date(:,2);
    simulp.year     = simulp.date(:,1);
    % simulp.year     = year(simul_date);
    

    % Dummy variable for week-days and week-end
    
    for i = 1 : n
        date = [num2str(simulp.year(i,1)),'-',num2str(simulp.month(i,1)),'-', num2str(simulp.day(i,1))];
        if weekday(date) == 1 || weekday(date) == 7
           simulp.weekend(i,1) = 1 ;
           simulp.weekday(i,1) = 0 ;
        else
           simulp.weekend(i,1) = 0 ;
           simulp.weekday(i,1) = 1 ;
        end
    end
   
    % Dummy variables for Peak / Off-Peaks hours
    
%     simulp.wd_offpk   = zeros(size(simulp.date,1),1) ;
%     simulp.wd_pk      = zeros(size(simulp.date,1),1) ;
%     
%     for i = 1 : size(simulp.hour,1) - 1;
%         if (simulp.hour(i,1) >= 7) && (simulp.hour(i,1) <= 19)  
%            simulp.wd_pk(i,1)    = simulp.weekday(i,1);
%         else
%            simulp.wd_offpk(i,1)	= simulp.weekday(i,1);
%         end
%     end


   simulp.prices1       = zeros(size(simulp.date,1),m);
   simulp.prices1(1,:)  = input2.prices (1,:);
% 
%         for i = 2 : size(simulp.prices1,1)
%             simulp.prices1(i,:) = simulp.weekday(i,1)   .* revrate_p_wd    .* (revlevel_p_wd    - simulp.prices1(i-1,:)) +...
%                                   simulp.weekend(i,1)   .* revrate_wend    .* (revlevel_wend - simulp.prices1(i-1,:)) +...
%                                   simulp.weekday(i,1)   .* sigma_p_wd      .* rnd(i,:) +...
%                                   simulp.weekend(i,1)   .* sigma_wend      .* rnd(i,:);  
%                              
%             simulp.prices1(i,:) = simulp.prices1(i,:) + simulp.prices1(i-1,:);
%         end
        for i = 2 : size(simulp.prices1,1)
            simulp.prices1(i,:) =  revrate_p_wd    .* ( revlevel_p_wd - simulp.prices1(i-1,:)) +...
                                   sigma_p_wd      .* rnd(i,:);
            simulp.prices1(i,:) = simulp.prices1(i,:) + simulp.prices1(i-1,:);
        end
        
   simulp.prices2   = zeros(size(simulp.date,1),m);
   simulp.prices2(1,:) = input2.prices (1,:);

        for i = 2 : size(simulp.prices2,1)-1
            simulp.prices2(i,:) =  revrate_p_wd    .* ( revlevel_p_wd    - simulp.prices2(i-1,:))+...
                                   sigma_p_wd      .* rnd(i,:)          + d_YFIT(i,:) ;
            simulp.prices2(i,:) = simulp.prices2(i,:) + simulp.prices2(i-1,:);
        end
        
        clear a_pk b_pk sig_pk m_pk gam_pk lamb_pk a_offpk b_offpk...
                sig_offpk m_offpk gam_offpk lamb_offpk a_weekend b_weekend...
                sig_weekend m_weekend gam_weekend lamb_weekend
            
       for i = 1 : 100 
            
%             [a_wend(i,1),b_wend(i,1),sig_wend(i,1),m_wend(i,1),gam_wend(i,1),lamb_wend(i,1)] =...
%                 mrjd_mle(nonzeros(simulp.weekend .*simulp.prices2(:,i)));
            
            [a_wd(i,1),b_wd(i,1),sig_wd(i,1),m_wd(i,1),gam_wd(i,1),lamb_wd(i,1)] =...
                     mrjd_mle(nonzeros(simulp.prices2(:,i)));
            
%             [a_weekend(i,1),b_weekend(i,1),sig_weekend(i,1),m_weekend(i,1),gam_weekend(i,1),lamb_weekend(i,1)] =...
%                 mrjd_mle(nonzeros(simulp.weekend .*simulp.prices2(:,i)));
       end
       
       % 
       
       disp('--------------------------------------------------------------------------------------------')
       disp('COMPARATIVE STUDY OF ACTUAL PARAMETERS VS. SIMULATION --------------------------------------')
       disp('--------------------------------------------------------------------------------------------')
       formatSpec1 = 'Stdev is in simulation wd: %1.4f  \n';
       formatSpec2 = 'Stdev is in actual wd: %1.4f \n';
       disp('--------------------------------------------------------------------------------------------')
       fprintf(formatSpec1,mean(sig_wd))
       fprintf(formatSpec2,sigma_p_wd)
       disp('--------------------------------------------------------------------------------------------')
       formatSpec1 = 'Intensity is in simulation wd: %1.4f \n';
       formatSpec2 = 'Intensity is in actual wd: %1.4f \n';
       disp('--------------------------------------------------------------------------------------------')
       fprintf(formatSpec1,mean(lamb_wd))
       fprintf(formatSpec2,lambda_p_wd)
       disp('--------------------------------------------------------------------------------------------')
       formatSpec1 = 'Mu (jump-size) is in simulation wd: %1.4f \n';
       formatSpec2 = 'Mu (jump-size) is in actual wd: %1.4f  \n';
       disp('--------------------------------------------------------------------------------------------')
       fprintf(formatSpec1,mean(m_wd))
       fprintf(formatSpec2,mu_p_wd)
       disp('--------------------------------------------------------------------------------------------')
       formatSpec1 = 'Gamma (jump-vol) is in simulation wd: %1.4f \n';
       formatSpec2 = 'Gamma (jump-vol) is in actual wd: %1.4f  \n';
       disp('--------------------------------------------------------------------------------------------')
       fprintf(formatSpec1,mean(gam_wd))
       fprintf(formatSpec2,gamma_p_wd)
       disp('--------------------------------------------------------------------------------------------')
       
       % Final simulation where we add the no-wind JUMP measure
       
       %numjumps_wend       = random('Poisson',max(0,lambda_wend-mean(lamb_wend)),m,1);
%        numjumps_wd    = random('Poisson',max(0,lambda_p_wd-mean(lamb_wd)),m,1);
% 
%        simulp.prices3   = simulp.prices2 ;
% 
%        for j = 1 : m
%             if numjumps_wd(j,1)>0 
%                 
%                 jump.place_wd    =  sort(unidrnd(size((simulp.prices2),1),     numjumps_wd(j,1),1));
%                 % jump.place_offpk =  sort(unidrnd(size((simulp.wd_offpk),1),  numjumps_offpk(j,1),1));
%                 % jump.place_wend  =  sort(unidrnd(size((simulp.weekend),1),   numjumps_wend(j,1),1));
%                 
%                 jump.size_wd  =  normrnd ( max(0,mu_p_wd-mean(m_wd)), max(0,gamma_p_wd-mean(gam_wd)), numjumps_wd(j,1), 1);
%                 % jump.size_offpk =  normrnd ( max(0,mu_offpk-mean(m_offpk)),   max(0,gamma_offpk-mean(gam_offpk)), numjumps_offpk(j,1), 1);
%                 % jump.size_wend  =  normrnd ( max(0,mu_wend-mean(m_weekend)),  max(0,gamma_wend-mean(gam_weekend)), numjumps_wend(j,1), 1);
%               
%                 for ii = 1 : size(jump.place_wd,1)
%                     while jump.size_wd(ii,1)>0 && simulp.prices2(jump.place_wd(ii,1))*jump.size_wd(ii,1) == 0
%                         jump.place_wd(ii,1) = jump.place_wd(ii,1) + 1;
%                     end
%                 end
%                 
%                 for ii = 1 : size(jump.place_offpk,1)
%                     while jump.size_offpk(ii,1)>0 && simulp.wd_offpk(jump.place_offpk(ii,1))*jump.size_offpk(ii,1) == 0
%                         jump.place_offpk(ii,1) = jump.place_offpk(ii,1) + 1;
%                     end
%                 end
%                 
% %                 for ii = 1 : size(jump.place_wend,1)
% %                     while jump.size_wend(ii,1)>0 && simulp.wd_wend(jump.place_wend(ii,1))*jump.size_wend(ii,1) == 0
% %                         jump.place_wend(ii,1) = jump.place_wend(ii,1) + 1;     
% %                     end
% %                 end
%                 
%                 jump_place = [jump.place_pk;    jump.place_offpk;   ];
%                 jump_size  = [jump.size_pk;     jump.size_offpk;    ];
%                 
%                 for i = 2 : size(simulp.prices3,1)
%                    for ii = 1 : size(jump_place,1);
%                       % if     i < jump_place(ii,1)
%                       %      simulp.prices3(i,j) = simulp.prices3(i,j) ;
%                       if i == jump_place(ii,1)
%                             simulp.prices3(i,j) = simulp.prices3(i,j) + jump_size(ii,1); 
%                       elseif i > jump_place(ii,1)
%                             simulp.prices3(i,j) = simulp.prices3(i-1,j) +...
%                                                   simulp.wd_pk(i,1)      .* revrate_pk       .* (revlevel_pk    - simulp.prices3(i-1,j))+...
%                                                   simulp.wd_offpk(i,1)   .* revrate_offpk    .* (revlevel_offpk - simulp.prices3(i-1,j))+...
%                                                                                                     simulp.wd_pk(i,1)      .* sigma_pk         .* rnd(i,j)+...
%                                                   simulp.wd_offpk(i,1)   .* sigma_offpk      .* rnd(i,j)+...
%                                                   (d_YFIT(i-1,j)) ; 
%                            
%                        end
%                    end 
%                 end
%             end
%        end

        
                                  
       
figure(9)

    subplot(4,1,1)
    plot(output.seasonal ,'-r','LineWidth',1); hold on
    plot(final_simul(1:size(d_YFIT,1),1),'-g','LineWidth',2); grid on; 
    xlabel('Hours in a year')
    ylabel('(Mwh)')
    title({'Simulation of hourly Wind Production (in green)'})

    subplot(4,1,2)
    plot(x1fit(1:size(d_YFIT,1),1),'-b','LineWidth',1); grid on; 
    xlabel('Hours in a year')
    ylabel('(Mwh)')
    title({'Simulation of hourly variation of Wind Production (in blue)'})
    
    subplot(4,1,3)
    plot(d_YFIT(:,1),'-r','LineWidth',2); grid on
    xlabel('Hours in year')
    ylabel('(GBP/MWH)')
    title({'Simulation of impact in Spot Prices produced by changes in Wind (in red)'})
    
    subplot(4,1,4)
    plot(simulp.prices2(:,1),'-r','LineWidth',2); grid on ; hold on;
    plot(simulp.prices1(:,1),'-b','LineWidth',1);hold on;
    plot(simulp.prices2(:,1)-simulp.prices1(:,1),'-k','LineWidth',1)
    xlabel('Hours in year')
    ylabel('(GBP/MWH)')
    title({'Simulation of Spot Prices produced by changes in Wind (in black)'})

figure(10)
   x = find(sum(simulp.prices1-simulp.prices2)>0);
   sim = 1;
   for i = 1 : size(x,2)
       subplot(3,5,sim)
%            plot(simulp.prices2(:,i),'-r','LineWidth',2); grid on ; hold on;
%            plot(simulp.prices3(:,i),'-k','LineWidth',1); hold on;
%            plot(simulp.prices1(:,i),'-b','LineWidth',1);
           plot(simulp.prices1(:,i)-simulp.prices2(:,i),'-g','LineWidth',1); hold on;
%            plot(simulp.prices2(:,i)-simulp.prices3(:,i),'-r','LineWidth',1); hold on;
           
           if sim < 15
              sim = sim + 1 ;
           elseif sim >= 15 ;    
              break ;
           end
   end
        

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
figure(11)
    subplot(3,1,1)
    hist(input2.prices)

    subplot(3,1,2)
    hist(simulp.prices1(:,1))
    
    subplot(3,1,3)
    hist(simulp.prices2(:,1))
    
    
% ---------------------------------------------------------------------------    
%% 9.  EXAMPLE of Pricing a Wind-Structured Product
% ---------------------------------------------------------------------------
    
 % Wind-Structured Product Calendar-Time
 
    firstdate   = '01-Jul-2014 ';
    lastdate    = '30-Jun-2028 ';
    simul_date  = [datenum(firstdate):datenum(lastdate)]';
    simulp.windfarm = datevec(simul_date,'yyyy');
    
    firstperiod_1  = '01-Jul-2014';
    firstperiod_2  = '30-Jun-2018';    
    simul_period1  = [datenum(firstperiod_1):datenum(firstperiod_2)]';
    simulp.period1 = datevec(simul_period1,'yyyy');

    secondperiod_1  = '01-Jul-2018';
    secondperiod_2  = '30-Jun-2023';
    simul_period2  = [datenum(secondperiod_1):datenum(secondperiod_2)]';
    simulp.period2 = datevec(simul_period2,'yyyy');
    
    thirdperiod_1  = '01-Jul-2023';
    thirdperiod_2  = '30-Jun-2028';
    simul_period3  = [datenum(thirdperiod_1):datenum(thirdperiod_2)]';
    simulp.period3 = datevec(simul_period3,'yyyy');
    
 % Volume Simulated for a wind-farm in Mwh
 
    maxprod  = 150.0 ; % Mwh/h <- as example
    simul    = 25 .* ( final_simul ./ maxprod ) .* maxprod ;
    
 % Payoff Prices per Mwh
    
    payoff_total1   = zeros(size(simulp.windfarm,1),m);
    payoff_cap1     = zeros(size(simulp.windfarm,1),m);
    payoff_floor1   = zeros(size(simulp.windfarm,1),m);
    
    payoff_total2   = zeros(size(simulp.windfarm,1),m);
    payoff_cap2     = zeros(size(simulp.windfarm,1),m);
    payoff_floor2   = zeros(size(simulp.windfarm,1),m);
    
    AIF            = zeros(size(simulp.windfarm,1),m);
    FIF            = zeros(size(simulp.windfarm,1),m);
    
    simulp.original1 = simulp.prices1;
    simulp.prices1   = simulp.original1(size(simulp.original1,1)-size(simulp.windfarm,1)+1:end,:);
    
    simulp.original2 = simulp.prices2;
    simulp.prices2   = simulp.original2(size(simulp.original2,1)-size(simulp.windfarm,1)+1:end,:);
    

    
    years = [(2014 : 2028);  (0.05 : 0.01 : 0.2)]';
    
    cap_strikes     = [55.55 56.66 57.79 58.95 60.13 61.33 62.56 63.81 65.08 66.38 69.57 72.83 74.29 75.78 77.29 78.84]'; 
    floor_strikes   = [ 30.30 30.90 31.52 32.15 30.04 27.88 28.43 29.00 29.58 30.17 27.67 25.11 25.62 26.13 26.65 27.19 ]';
    
    for j = 1 : m 
        
        for i = 1 : size(simulp.period1,1)
            
            if      (simulp.windfarm(i,1) == simulp.period1(i,1)) && ...
                    (simulp.windfarm(i,2) == simulp.period1(i,2)) && ...
                    (simulp.windfarm(i,3) == simulp.period1(i,3))

                    [x,y] = find (simulp.windfarm(i,1) == years);

                    payoff_cap1(i,j)          = max(0 , simulp.prices1(i,j) - cap_strikes(x,1) ).* 0.9;
                    payoff_floor1(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices1(i,j) ).* 0.9;

                    payoff_cap2(i,j)          = max(0 , simulp.prices2(i,j) - cap_strikes(x,1) ).* 0.9;
                    payoff_floor2(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices2(i,j) ).* 0.9;                    
                    
                    AIF(i,j)    = years(x,2) .* simulp.prices1(i,j);
                    FIF(i,j)    = simulp.prices1(i,j).* 0.9;
                    
                    payoff_total1(i,j)   = simulp.prices1(i,j).* 0.9 - payoff_cap1(i,j) +  payoff_floor1(i,j) ...
                                            - years(x,2) .* simulp.prices1(i,j) ;
                                        
                    payoff_total2(i,j)   = simulp.prices2(i,j).* 0.9 - payoff_cap2(i,j) +  payoff_floor2(i,j) ...
                                            - years(x,2) .* simulp.prices2(i,j) ;
            end
            
        end
        
        for i = size(simulp.period1,1)+1 : size(simulp.period1,1)+size(simulp.period2,1)
            if      (simulp.windfarm(i,1) == simulp.period2(i-size(simulp.period2,1),1)) && ...
                    (simulp.windfarm(i,2) == simulp.period2(i-size(simulp.period2,1),2))&& ...
                    (simulp.windfarm(i,3) == simulp.period2(i-size(simulp.period2,1),3))                                   

                    [x,y] = find (simulp.windfarm(i,1) == years);
                    
                    payoff_cap1(i,j)          = max(0 , simulp.prices1(i,j) - cap_strikes(x,1) ).* 0.8;
                    payoff_floor1(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices1(i,j) ).* 0.8;

                    AIF(i,j) = years(x,2) .* simulp.prices1(i,j);
                    FIF(i,j) = simulp.prices1(i,j).* 0.8;
                    
                    payoff_cap2(i,j)          = max(0 , simulp.prices2(i,j) - cap_strikes(x,1) ).* 0.8;
                    payoff_floor2(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices2(i,j) ).* 0.8; 
                    
                    payoff_total1(i,j)   = simulp.prices1(i,j).* 0.8 - payoff_cap1(i,j) +  payoff_floor1(i,j) ...
                                            - years(x,2) .* simulp.prices1(i,j) ;
                                        
                    payoff_total2(i,j)   = simulp.prices2(i,j).* 0.8 - payoff_cap2(i,j) +  payoff_floor2(i,j) ...
                                            - years(x,2) .* simulp.prices2(i,j) ;
            end
        end
        
        for i = size(simulp.period1,1)+size(simulp.period2,1)+1 : size(simulp.period1,1)+size(simulp.period2,1)+size(simulp.period3,1)-1
            if      (simulp.windfarm(i,1) == simulp.period3(i-size(simulp.period3,1)-size(simulp.period2,1)+1,1)) && ...
                    (simulp.windfarm(i,2) == simulp.period3(i-size(simulp.period3,1)-size(simulp.period2,1)+1,2)) && ...
                    (simulp.windfarm(i,3) == simulp.period3(i-size(simulp.period3,1)-size(simulp.period2,1)+1,3))
                                    
                    [x,y] = find (simulp.windfarm(i,1) == years);
                    
                    payoff_cap1(i,j)          = max(0 , simulp.prices1(i,j) - cap_strikes(x,1) ).* 0.7;
                    payoff_floor1(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices1(i,j) ).* 0.7;

                    AIF(i,j) = years(x,2) .* simulp.prices1(i,j);
                    FIF(i,j) = simulp.prices1(i,j).* 0.7;
                    
                    payoff_cap2(i,j)          = max(0 , simulp.prices2(i,j) - cap_strikes(x,1) ).* 0.7;
                    payoff_floor2(i,j)        = max(0 , floor_strikes(x,1) - simulp.prices2(i,j) ).* 0.7; 
                    
                    payoff_total1(i,j)   = simulp.prices1(i,j).* 0.7 - payoff_cap1(i,j) +  payoff_floor1(i,j) ...
                                            - years(x,2) .* simulp.prices1(i,j) ;
                                        
                    payoff_total2(i,j)   = simulp.prices2(i,j).* 0.7 - payoff_cap2(i,j) +  payoff_floor2(i,j) ...
                                            - years(x,2) .* simulp.prices2(i,j) ;
            end
        end
            
    end
    
    
    % Case 1: Wind Production does not affect prices
    
    final_value_per_Mwh1            = sum (mean(payoff_total1,2));
    final_value_pMwh_prctile1       = sum (prctile(payoff_total1,[5 50 95],2));

    production_vector               = [ 58.3 58.3 58.3 38.8 38.8 38.8 38.8 38.8 38.8 58.3 58.3 58.3 ];
    
    proxy_old_valuation_low1        = sum (production_vector .* 1000 .* mean (prctile (payoff_total1,[5],2)));
    proxy_old_valuation_mean1       = sum (production_vector .* 1000 .* mean (prctile (payoff_total1,[50],2)));
    proxy_old_valuation_high1       = sum (production_vector .* 1000 .* mean (prctile (payoff_total1,[95],2)));
    
    final_simulation_value1         = payoff_total1 .* London_simul( size(London_simul,1) - size(payoff_total1) + 1 : end,: );
    final_simulation_percentiles1   = sum(prctile(final_simulation_value1(1:end-1,:),[5 50 95],2));

    new_value_wind1                 = sum(mean(final_simulation_value1,2));
    
    
    % Case 2: Wind Production affects prices
    
    final_value_per_Mwh2            = sum(mean(payoff_total2,2));
    final_value_pMwh_prctile2       = sum(prctile(payoff_total2,[5 50 95],2));

    proxy_old_valuation_low2        = sum (production_vector .* 1000 .* mean(prctile(payoff_total2,[5],2)));
    proxy_old_valuation_mean2       = sum (production_vector .* 1000 .* mean(prctile(payoff_total2,[50],2)));
    proxy_old_valuation_high2       = sum (production_vector .* 1000 .* mean(prctile(payoff_total2,[95],2)));
    
    final_simulation_value2         = payoff_total2 .* London_simul(size(London_simul,1) - size(payoff_total2) +1 : end,:);
    final_simulation_percentiles2   = sum (prctile (final_simulation_value2 (1:end-1,:),[5 50 95],2));
    
    new_value_wind2                 = sum(mean(final_simulation_value2,2));
    

    
    % Screen Report

    disp('------------------------------------------------------------------------------------------------------')
    disp (' 9. PRICING Windfarm STRUCTURE ' )
    disp('------------------------------------------------------------------------------------------------------')
    disp('  DETERMINISTIC WIND VALUATION ')
    disp (' ')
    disp('  -> windfarm value 1 Mwh (15 years contract): ')    
    disp (' ')         
    sentence1 = '\t %10.1f GBP  with 95pc of interval ( %10.1f , %10.1f ) \n';
    fprintf(sentence1, final_value_per_Mwh1, final_value_pMwh_prctile1(1,1), final_value_pMwh_prctile1(1,3))
    disp (' ')
    disp('  -> windfarm value 500.Gwh (15 years contract) according with calendar production in the Term Sheet: ')
    disp (' ')
    sentence2 = '\t %10.1f GBP  with 95pc of interval ( %10.1f , %10.1f ) \n';
    fprintf(sentence2, proxy_old_valuation_mean1, proxy_old_valuation_low1, proxy_old_valuation_high1)
    disp(' ')
    disp('------------------------------------------------------------------------------------------------------')
    disp('  STOCHASTIC WIND VALUATION ')
    disp(' ')
    disp (' CASE 1 : Wind production does NOT affect prices ')
    disp (' ')
    disp('  -> windfarm value with Stochastic Wind Production of the Wind Farm:  ');
    disp (' ')
    sentence3 = '\t %10.1f GBP  with 95pc of interval ( %10.1f , %10.1f ) \n';
    fprintf(sentence3, new_value_wind1, final_simulation_percentiles1(1,1), final_simulation_percentiles1(1,3))
    disp (' ')
    disp (' ')
    disp (' CASE 2 : Wind production AFFECTS prices using fundamental model')
    disp (' ')
    disp('  -> windfarm value with Stochastic Wind Production of the Wind Farm:  ');
    disp(' ')
    fprintf(sentence3, new_value_wind2, final_simulation_percentiles2(1,1), final_simulation_percentiles2(1,3))
    disp (' ')
    disp (' ')
    disp('------------------------------------------------------------------------------------------------------')
    
    
    
    
    
    
    
    figure(14)
    
    subplot(2,2,1)
    plot(payoff_total1)
    % plot(mean(payoff_total,2))    
    title('Simulations of Payoff total for 1Mwh simulated')
    grid on
    
    subplot(2,2,2)
    
    plot(prctile(payoff_total1(1:end-1,:),[5 50 95],2))
    title('Percentiles of Daily final Payoff for 1Mwh')
    grid on
    
    subplot(2,2,3)
    plot(mean(payoff_cap1(1:end-1,:),2),'-g')
    title('Mean of Daily Cap Value for 1Mwh')
    grid on
    
    subplot(2,2,4)
    plot(mean(payoff_floor1(1:end-1,:),2),'-r')
    title('Mean of Daily Floor Value for 1Mwh')
    grid on
    

    figure(15)
    
    subplot(1,2,1)
    plot(final_simulation_value1(1:end-1,:))
    title('Stochastic Production * stochastic Payoff per Mwh')
    grid on
    
    subplot(1,2,2)
%     plot(mean(final_simulation_value,2))
%     hold on
    plot(prctile(final_simulation_value1(1:end-1,:),[5 50 95],2))
    title('Percentiles 5-50-95% of payoffs')
    grid on





    
    