function CV = one_binding_site_model(p, L1, L2, ligandchoice)
    % Inputs:   p - parameter vector
    %           L1 : background concentration of Lasp
    %           L2 : background concentration of meAsp
    %           ligandchoice : 0 output Lasp khalf CV, 1 output meAsp khalf
    %           CV
    
    % Parameters
    Ki1 = p(1); % Dissociation constant for Lasp at site 1
    Ki2 = p(2); % Dissociation constant for meAsp at site 1
    K0_meAsp =  p(3); % Mean Khalf for meAsp in 0-background
    sigma0_meAsp = p(4); % std of Khalf for meAsp in 0-background
    K0_Lasp = p(5); % Mean Khalf for Lasp in 0-background
%     K0_Lasp = 0.0799; % Mean Khalf for Lasp in 0-background
    sigma0_lAsp = p(6); % standard deviation L-asp Khalf in 0-background
%     CV0 = p(4);
%     c = p(7);

%     K0_meAsp =  1.2552; % Mean Khalf for meAsp in 0-background
%     Ki2 = 18;
    
    % Calculate Khalf CVs
    if ligandchoice == 0
        CV = sigma0_lAsp./(K0_Lasp + L1./(1+L1./Ki1 + L2./Ki2));
%         CV = CV0.*K0_Lasp./(K0_Lasp + L1./(1+L1./Ki1 + L2./Ki2));
    elseif ligandchoice == 1
        CV = sigma0_meAsp./(K0_meAsp + L2./(1+L1./Ki1 + L2./Ki2));
%         CV = CV0.*K0_meAsp./(K0_meAsp + L2./(1+L1./Ki1 + L2./Ki2));
    end
    return
end