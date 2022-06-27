function settings = synths(settings)

settings.A0 = [];
settings.Af = [];
settings.V = [];
settings.L = [];
settings.r = {};
settings.x0 = [];
settings.xf = [];

%% PITOT
for i = 1:length(settings.Lpitot)
    settings.A0(end+1) = pi * settings.Dpitot(i)^2/4;
    settings.Af(end+1) = pi * settings.Dpitot(i+1)^2/4;
    settings.L(end+1) = settings.Lpitot(i);
    
    r = @(x) (settings.Dpitot(i+1)-settings.Dpitot(i))/(2*settings.Lpitot(i)) * x + settings.Dpitot(i)/2;
    f = @(x) r(x).^2;
    settings.V(end+1) = pi * integral(f, 0, settings.Lpitot(i), 'AbsTol', 1e-6, 'RelTol', 1e-6);
    settings.r{end+1} = @(x) r(x);
    settings.x0(end+1) = 0;
    settings.xf(end+1) = settings.Lpitot(i);
end

%% NOSECONE
if settings.Dnose(1) == 0
    settings.A0(end+1) = 0;
else
    settings.A0(end+1) = pi * settings.Dnose(1)^2/4;
end
settings.Af(end+1) = pi * settings.Dnose(2)^2 / 4;

switch settings.Tnose
    case 'CONIC'
        r = @(x) settings.Dnose(2)/(2*settings.LnoseTrue) * x;
        xx = fzero(@(x) r(x) - settings.Dnose(1)/2, 0.01);
        if settings.Dnose(1)/2 == 0
            xx = 0;
            settings.Lnose = settings.LnoseTrue;
        else
            xx = fzero(@(x) r(x) - settings.Dnose(1)/2, 0.01);
            settings.Lnose = settings.LnoseTrue - xx;
        end
        f = @(x) (r(x) < settings.Dnose(1)/2).*0 + (r(x) >= settings.Dnose(1)/2).*r(x).^2;
        settings.V(end+1) = pi * integral(f, 0, settings.LnoseTrue, 'AbsTol', 1e-6, 'RelTol', 1e-6);
        settings.r{end+1} = @(x) r(x + xx);
        settings.x0(end+1) = 0;
        settings.xf(end+1) = settings.Lnose;
        settings.xx = xx;
        
    case 'HAACK'
        theta = @(x) acos(1 - 2*x/settings.LnoseTrue);
        rr = @(x) settings.Dnose(2)/(2*sqrt(pi)) * sqrt(theta(x) - sin(2*theta(x))/2 + (0)*sin(theta(x)).^3);
        r = @(x) (x<0).*0 + (x>=0).*rr(x);
        if settings.Dnose(1)/2 == 0
            xx = 0;
            settings.Lnose = settings.LnoseTrue;
        else
            xx = fzero(@(x) r(x) - settings.Dnose(1)/2, 0.01);
            settings.Lnose = settings.LnoseTrue - xx;
        end
        f = @(x) (r(x) < settings.Dnose(1)/2).*0 + (r(x) >= settings.Dnose(1)/2).*r(x).^2;
        f1 = @(x) r(x+xx).^2;
        settings.V(end+1) = pi * integral(f1, 0, settings.Lnose, 'AbsTol', 1e-6, 'RelTol', 1e-6);
        settings.r{end+1} = @(x) r(x+xx) ;
        settings.x0(end+1) = 0;
        settings.xf(end+1) = settings.Lnose;
        settings.xx = xx;
end
settings.L(end+1) = settings.Lnose;

%% CENTERBODY
settings.A0(end+1) = pi * settings.Dcentr^2 / 4;
settings.Af(end+1) = pi * settings.Dcentr^2 / 4;
settings.L(end+1)  = settings.Lcentr;
settings.V(end+1)  = settings.A0(end) * settings.Lcentr;
settings.r{end+1} = @(x) settings.Dcentr/2 .* ones(size(x));
settings.x0(end+1) = 0;
settings.xf(end+1) = settings.Lcentr;

%% AFTERBODY
settings.A0(end+1) = pi * settings.Daft(1)^2 / 4;
settings.Af(end+1) = pi * settings.Daft(2)^2 / 4;
settings.L(end+1)  = settings.Laft;
r = @(x) (settings.Daft(2)-settings.Daft(1))/(2*settings.Laft) * x + settings.Daft(1)/2;
f = @(x) r(x).^2;
settings.V(end+1) = pi * integral(f, 0, settings.Laft, 'AbsTol', 1e-6, 'RelTol', 1e-6);
settings.r{end+1} = @(x) r(x);
settings.x0(end+1) = 0;
settings.xf(end+1) = settings.Laft;

%% FINS
settings.Yfin = @(x) interp1(settings.xi, settings.yi, x);
settings.Cfin = @(y) interp1(settings.yi(1:floor(length(settings.yi)/2)), settings.xi(end:-1:floor(length(settings.yi)/2)+1) - settings.xi(1:floor(length(settings.yi)/2)), y);
settings.YLEfin = @(y) interp1(settings.yi(1:floor(length(settings.yi)/2)), settings.xi(1:floor(length(settings.yi)/2)), y);
settings.xf(end+1) = max(settings.xi);

Afin = integral(settings.Yfin, 0, max(settings.xi), 'AbsTol', 1e-6, 'RelTol', 1e-6);
settings.Afin = Afin;
f = @(x) settings.Cfin(x).^2;
settings.c_bar = 1/Afin * integral(f, 0, max(settings.yi), 'AbsTol', 1e-6, 'RelTol', 1e-6);
f = @(x) settings.YLEfin(x) .* settings.Cfin(x);
settings.yMAC = 1/Afin * integral(f, 0, max(settings.yi), 'AbsTol', 1e-6, 'RelTol', 1e-6);
settings.s = max(settings.yi);
settings.Tc = atan(((settings.xi(find(settings.yi == max(settings.yi), 1, 'first'))) + (settings.Cfin(max(settings.yi)))/2 - settings.xi(end)/2) / max(settings.yi));


