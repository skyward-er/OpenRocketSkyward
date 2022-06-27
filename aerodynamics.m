function out = aerodynamics(settings, alpha)

Aref = settings.Aref;
d = settings.Dcentr;
Na = length(alpha);
Nf = settings.Nfins;

body.CN_alpha = zeros(Na, 2);
body.CN = zeros(Na, 2);
body.CM = zeros(Na, 2);
body.Cm_alpha = zeros(Na, 2);
body.CN_alpha0 = zeros(Na, 2);
body.CN_alphaL = zeros(Na, 2);
body.Cm = zeros(Na, 2);
body.X = zeros(Na, 2);
body.X0 = zeros(Na, 2);
body.XL = zeros(Na, 2);

%% BODY
for k = 1:Na
    for i = 1:length(settings.A0)
        Af = settings.Af(i);
        A0 = settings.A0(i);
        l = settings.L(i);
        V = settings.V(i);

        %% CN
        if alpha(k) == 0
            body.CN_alpha0(k, i) = (2/Aref) * (Af - A0);
        else
            body.CN_alpha0(k, i) = (2/Aref) * (Af - A0) * (sin(alpha(k))/alpha(k));
        end

        %% CM
        if alpha(k) == 0
            body.Cm_alpha(k, i) = (2/(Aref*d)) * (l * Af - V);
        else
            body.Cm_alpha(k, i) = (2/(Aref*d)) * (l * Af - V) * (sin(alpha(k))/alpha(k));
        end
        body.Cm(k, i) = body.Cm_alpha(k, i) * alpha(k);

        %% X
        f = @(x) x.*(2*settings.r{i}(x));
        f1 = @(x) settings.r{i}(x);

        if alpha(k) == 0
            body.CN_alphaL(k, i) = 1.1 * 2*(2*integral(f1, 0, settings.xf(i), 'AbsTol', 1e-6, 'RelTol', 1e-6)/Aref)*sin(alpha(k))*cos(alpha(k));
        else
            body.CN_alphaL(k, i) = 1.1 * (2*integral(f1, 0, settings.xf(i), 'AbsTol', 1e-6, 'RelTol', 1e-6)/Aref)*sin(alpha(k))^2/alpha(k);
        end
        body.CN_alpha(k, i) = body.CN_alpha0(k, i) + body.CN_alphaL(k, i);
            
        body.CN(k, i) = alpha(k) * body.CN_alpha(k, i);
        body.CM(k, i) = alpha(k) * body.Cm_alpha(k, i);
        
        body.XL(k, i) = integral(f, 0, settings.xf(i), 'AbsTol', 1e-6, 'RelTol', 1e-6) / (2*integral(f1, 0, settings.xf(i), 'AbsTol', 1e-6, 'RelTol', 1e-6)) + sum(settings.L(1:i-1));
        
        if settings.r{i}(0) == settings.r{i}(settings.xf(i))
            body.X(k, i) = body.XL(k, i);
        else
            body.X0(k, i) = body.Cm_alpha(k, i)/body.CN_alpha0(k, i) * d  + sum(settings.L(1:i-1));
            body.X(k, i) = (body.X0(k, i)*body.CN_alpha0(k, i) + ...
                body.CN_alphaL(k, i)*body.XL(k, i))/(body.CN_alpha0(k, i) + body.CN_alphaL(k, i));
        end
    end
    
    body.CP(k) = sum(body.X(k, :).*body.CN_alpha(k, :))/sum(body.CN_alpha(k, :));
end

out.body = body;

%% FINS
fins.CN_alphai = zeros(Na, Nf);
fins.CN_alpha = zeros(Na, 1);

for k = 1:Na
    for i = 1:Nf
        beta = sqrt(1 - settings.Mach^2);
        fins.CN_alphai(k, i) = (2*pi*settings.s^2/Aref)/...
            (1 + sqrt(1 + ((beta*settings.s^2)/(settings.Afin*cos(settings.Tc)))^2));
    end
    tau = 1 + settings.Dcentr/2 /(settings.s + settings.Dcentr/2);
    fins.CN_alpha(k) = fins.CN_alphai(k, 1) * Nf/2 * 1 * tau;
    fins.X = sum(settings.Lpitot) + settings.Lnose + settings.Xle + settings.yMAC + 0.25*settings.c_bar;
end

out.fins = fins;

%% TOTAL
total.CN_alpha = zeros(Na, 1);
total.X = zeros(Na, 1);

for i = 1:Na
    for k1 = 1:size(body.CN_alpha, 2)
        X_body = sum(body.X(i, :).*body.CN_alpha(i, :))/sum(body.CN_alpha(i, :));
        CN_alpha_body = sum(body.CN_alpha(i, :));
    end
    
    for k2 = 1:size(fins.CN_alpha, 2)
        total.X(i) = (X_body*CN_alpha_body + fins.X*fins.CN_alpha(i, k2))/(CN_alpha_body + fins.CN_alpha(i, k2));
        total.CN_alpha(i) = CN_alpha_body + fins.CN_alpha(i, k2);
    end
end


out.total = total;


%% OUT
% out.CN_alpha = CN_alpha;
% out.CN_alpha0 = CN_alpha0;
% out.CN_alphaL = CN_alphaL;
% out.CN = CN;
% out.CM = CM;
% out.Cm_alpha = Cm_alpha;
% out.CP = out.CP';
% out.X = X;
% out.XL = XL;
% out.X0 = X0;
% out.alpha = alpha;

end