function CV = one_binding_site_model(p, L1, L2, ligandchoice)
    % Inputs:   p - parameter vector
    %           L1 : background concentration of Lasp
    %           L2 : background concentration of meAsp
    %           ligandchoice : 0 output Lasp khalf CV, 1 output meAsp khalf
    %           CV
    
    % Parameters
    Ki1 = p(1); % Dissociation constant for Lasp at site 1
    Ki2 = p(2); % Dissociation constant for meAsp at site 1
    mu = p(3); % Mean of phenotype P0 distribution
    sigma = p(4); % standard deviation of phenotype distribution

    % Functions
    f = @(L1, L2) (1 + L1./Ki1 + L2./Ki2);

    % Khalf functions
    Khalf_Lasp = @(L1_0, L2_0, P0) (P0.*f(L1_0, L2_0) - L2_0./Ki2 - 1).*Ki1;
    Khalf_meAsp = @(L1_0, L2_0, P0) (P0.*f(L1_0, L2_0) - L1_0./Ki1 - 1).*Ki2;


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