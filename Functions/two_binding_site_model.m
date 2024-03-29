function CV = two_binding_site_model(p, L1, L2, ligandchoice)
    % Inputs:   p - parameter vector
    %           L1 : background concentration of Lasp
    %           L2 : background concentration of meAsp
    %           ligandchoice : 0 output Lasp khalf CV, 1 output meAsp khalf
    %           CV
    
    % Parameters
    Ki1 = p(1); % Dissociation constant for Lasp at site 1
    Ki2 = p(2); % Dissociation constant for meAsp at site 1
    Ki11 = p(3); % Dissociation constant for Lasp at site 2
    mu = p(4); % Mean of phenotype P0 distribution
    sigma = p(5); % standard deviation of phenotype distribution

    % Functions
    f = @(L1, L2) (1+L1./Ki1+L2./Ki2).*(1+L1./Ki11);
    A1 = 1./(Ki11.*Ki1);
    B1 = @(L2) 1./Ki11 + 1./Ki1 + L2./(Ki11.*Ki2);
    C1 = @(L1, L2, P0) 1 + L2./Ki2 - P0.*f(L1, L2);

    % Khalf functions
    Khalf_Lasp = @(L1_0, L2_0, P0) (-B1(L2_0) + sqrt(B1(L2_0).^2 - 4.*A1.*C1(L1_0, L2_0, P0)))./(2.*A1);
    Khalf_meAsp = @(L1_0, L2_0, P0) (P0.*f(L1_0, L2_0) ./ (1 + L1_0./Ki11) - (1 + L1./Ki1)).*Ki2;

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