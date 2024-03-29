function CV = negative_cooperativity_model(p, L1, L2, ligandchoice)
    % Inputs:   p - parameter vector
    %           L1 : background concentration of Lasp
    %           L2 : background concentration of meAsp
    %           ligandchoice : 0 output Lasp khalf CV, 1 output meAsp khalf
    %           CV
    
    % Parameters
    Ki1 = p(1);
    Ki2 = p(2);
    e11 = p(3);
    e12 = p(4);
    e22 = p(5);
    mu = p(6);
    sigma = p(7);
%     mu = 1.33;
%     sigma = 0.39;

    % Functions
    f = @(L1, L2) 1 + 2.*L1./Ki1 + 2.*L2./Ki2 + exp(-e11).*(L1./Ki1).^2 +...
        exp(-e22).*(L2./Ki2).^2 + exp(-e12).*(L1./Ki1).*(L2./Ki2);
    A1 = exp(-e11)./(Ki1.^2);
    B1 = @(L2) 2./Ki1 + exp(-e12).*(1./Ki1).*(L2./Ki2);
    C1 = @(L1, L2, P0) 1 + 2.*(L2./Ki2) + exp(-e22).*(L2./Ki2).^2 - P0.*f(L1, L2);
    
    A2 = exp(-e22)./(Ki2.^2);
    B2 = @(L1) 2./Ki2 + exp(-e12).*(1./Ki2).*(L1./Ki1);
    C2 = @(L1, L2, P0) 1 + 2.*(L1./Ki1) + exp(-e11).*(L1./Ki1).^2 - P0.*f(L1, L2);

    % Khalf functions
    Khalf_Lasp = @(L1_0, L2_0, P0) (-B1(L2_0) + sqrt(B1(L2_0).^2 - 4.*A1.*C1(L1_0, L2_0, P0)))./(2.*A1);
    Khalf_meAsp = @(L1_0, L2_0, P0) (-B2(L1_0) + sqrt(B2(L1_0).^2 - 4.*A2.*C2(L1_0, L2_0, P0)))./(2.*A2);
%     Khalf_meAsp = @(L1_0, L2_0, P0) -C2(L1_0, L2_0, P0)./B2(L1_0);

    % Sample phenotypes
    nSamples = 10000;
    rng(1)
    P0_samples = 1+lognrnd(mu, sigma, [nSamples, 1]);
    
    % Calculate Khalf CVs
    if ligandchoice == 0
        Khalf_temp = Khalf_Lasp(L1, L2, P0_samples);
        CV = std(Khalf_temp)./mean(Khalf_temp);
    elseif ligandchoice == 1
        Khalf_temp = Khalf_meAsp(L1, L2, P0_samples);
        CV = std(Khalf_temp)./mean(Khalf_temp);
    end
    return
end