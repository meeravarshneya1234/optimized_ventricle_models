% Ohara et al model
% 2011
% human ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.celltype = 'endo';
p.Nao=140.0;
p.Cao=1.8;
p.Ko=5.4;

%physical constants
p.R=8314.0;
p.T=310.0;
p.F=96485.0;
p.Cm=1.0;                     %uF

%cell geometry
p.L=0.01;
p.rad=0.0011;
p.vcell=1000*3.14*p.rad*p.rad*p.L;
p.Ageo=2*3.14*p.rad*p.rad+2*3.14*p.rad*p.L;
p.Acap=2*p.Ageo;
p.vmyo=0.68*p.vcell;
p.vnsr=0.0552*p.vcell;
p.vjsr=0.0048*p.vcell;
p.vss=0.02*p.vcell;

%jsr constants
p.bt=4.75;
p.a_rel=0.5*p.bt;

% computed quantities that do not change during simulation
c.GNa=75;
c.Gto=0.02;
c.GK1=0.1908;
c.GKb=0.003;
c.GpCa=0.0005;

c.PNab=3.75e-10;
c.PCab=2.5e-8;

c.SERCA_total = 1 ;
c.RyR_total = 1 ;

% % Parameters changed in the optimized model
c.GKs_=0.0034*5.75;
c.GKr_=0.046*1.00;
c.PCa_=0.0001*2.01;
c.Gncx=0.0008*2.95;
c.Pnak=30*9.12;
c.GNaL=0.0075*1.00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCL = 1000 ;  % Interval bewteen stimuli,[ms]
stim_delay = 100 ; % Time the first stimulus, [ms]
stim_dur = 2 ; % Stimulus duration
stim_amp = 20 ; % Stimulus amplitude 
nBeats = 100 ; % Number of beats to simulate 

stim_starts = stim_delay + PCL*(0:nBeats-1)  ;
stim_ends = stim_starts + stim_dur ;

% Create intervals for each beat 
simints = 3*nBeats ;
for i=1:nBeats
    intervals(3*i-2,:) = [PCL*(i-1),stim_starts(i)] ; %beginning 
    intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ; %stimulus 
    intervals(3*i,:) = [stim_ends(i),PCL*i] ; % stimulus ends 
end
tend = nBeats*PCL ;              % end of simulation, ms
intervals(end,:) = [stim_ends(end),tend] ;

% Determine when to apply stim_amp or 0 amp  
Istim = zeros(simints,1) ;
stimindices = 3*(1:nBeats) - 1 ; % apply stimulus on second part of intervals
Istim(stimindices) = -stim_amp ; 

ssbefore = 1;% Run the model with no stimulus for 60 seconds. 
numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ic.V=-87;
ic.Nai=7;
ic.Nass=ic.Nai;
ic.Ki=145;
ic.Kss=ic.Ki;
ic.Cai=1.0e-4;
ic.Cass=ic.Cai;
ic.Cansr=1.2;
ic.Cajsr=ic.Cansr;
ic.m=0;
ic.hf=1;
ic.hs=1;
ic.j=1;
ic.hsp=1;
ic.jp=1;
ic.mL=0;
ic.hL=1;
ic.hLp=1;
ic.a=0;
ic.iF=1;
ic.iS=1;
ic.ap=0;
ic.iFp=1;
ic.iSp=1;
ic.d=0;
ic.ff=1;
ic.fs=1;
ic.fcaf=1;
ic.fcas=1;
ic.jca=1;
ic.nca=0;
ic.ffp=1;
ic.fcafp=1;
ic.xrf=0;
ic.xrs=0;
ic.xs1=0;
ic.xs2=0;
ic.xk1=1;
ic.Jrelnp=0;
ic.Jrelp=0;
ic.CaMKt=0;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Ohara_opt;

statevar_i = y0;
%%% let model rest for 60 seconds
if (ssbefore)
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [~,posstatevars] = ode15s(odefcn,[0,60000],statevar_i,options,0,p,c) ;
    statevar_i = posstatevars(end,:) ;
end

%%% stimulate cell
if (nBeats > numbertokeep)
    for i=1:simints-3*numbertokeep
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        statevar_i = posstatevars(end,:) ;
        t = post(end) ;
    end % for
    statevars = statevar_i ;
    for i=simints-3*numbertokeep+1:simints
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end % for
else
    t = 0 ;
    statevars = statevar_i ;
    for i=1:simints
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end % for
end % if

t = t - min(t) ;
V = statevars(:,V_ind);
state_variables = num2cell(statevars,1) ;
APD= find_APD(t,V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 5:  Plot AP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t,V,'linewidth',2)
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('time (ms)')
ylabel('Voltage (mV)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 6:  Plot IKs, IKr, ICaL, INaL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V, Nai, Nass, Ki, Kss, Cai, Cass, Cansr, Cajsr, m, hf, hs, j, hsp, jp, mL, hL,...
    hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, ffp, fcafp,...
    xrf, xrs, xs1, xs2, xk1, Jrelnp, Jrelp, CaMKt] = deal(state_variables{:});

% % % ICaL
KmCaMK=0.15;
CaMKo=0.05;
KmCaM=0.0015;
CaMKb=CaMKo*(1.0-CaMKt)./(1.0+KmCaM./Cass);
CaMKa=CaMKb+CaMKt;
vffrt=V*p.F*p.F/(p.R*p.T);
vfrt=V*p.F/(p.R*p.T);
Aff=0.6;
Afs=1.0-Aff;
f=Aff*ff+Afs*fs;
Afcaf=0.3+0.6./(1.0+exp((V-10.0)/10.0));
Afcas=1.0-Afcaf;
fca=Afcaf.*fcaf+Afcas.*fcas;
fp=Aff.*ffp+Afs.*fs;
fcap=Afcaf.*fcafp+Afcas.*fcas;
PhiCaL=4.0*vffrt.*(Cass.*exp(2.0*vfrt)-0.341*p.Cao)./(exp(2.0*vfrt)-1.0);
PhiCaNa=1.0*vffrt.*(0.75*Nass.*exp(1.0*vfrt)-0.75*p.Nao)./(exp(1.0*vfrt)-1.0);
PhiCaK=1.0*vffrt.*(0.75*Kss.*exp(1.0*vfrt)-0.75*p.Ko)./(exp(1.0*vfrt)-1.0);
PCap=1.1*c.PCa_;
PCaNa=0.00125*c.PCa_;
PCaK=3.574e-4*c.PCa_;
PCaNap=0.00125*PCap;
PCaKp=3.574e-4*PCap;
fICaLp=(1.0./(1.0+KmCaMK./CaMKa));
ICa_L=(1.0-fICaLp).*c.PCa_.*PhiCaL.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCap.*PhiCaL.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaNa=(1.0-fICaLp).*PCaNa.*PhiCaNa.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaNap.*PhiCaNa.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaK=(1.0-fICaLp).*PCaK.*PhiCaK.*d.*(f.*(1.0-nca)+jca.*fca.*nca)+fICaLp*PCaKp.*PhiCaK.*d.*(fp.*(1.0-nca)+jca.*fcap.*nca);
ICaL = ICa_L + ICaNa + ICaK;

% % % IKr
EK=(p.R*p.T/p.F)*log(p.Ko./Ki);
Axrf=1.0./(1.0+exp((V+54.81)/38.21));
Axrs=1.0-Axrf;
xr=Axrf.*xrf+Axrs.*xrs;
rkr=1.0./(1.0+exp((V+55.0)/75.0))*1.0./(1.0+exp((V-10.0)/30.0));
IKr=c.GKr_*sqrt(p.Ko/5.4).*xr.*rkr.*(V-EK);

% % % IKs
PKNa=0.01833;
EKs=(p.R*p.T/p.F)*log((p.Ko+PKNa*p.Nao)./(Ki+PKNa.*Nai));
KsCa=1.0+0.6./(1.0+(3.8e-5./Cai).^1.4);
IKs=c.GKs_.*KsCa.*xs1.*xs2.*(V-EKs);

% % % INaL 
ENa=(p.R*p.T/p.F)*log(p.Nao./Nai);
CaMKo=0.05;
KmCaM=0.0015;
CaMKb=CaMKo.*(1.0-CaMKt)./(1.0+KmCaM./Cass);
CaMKa=CaMKb+CaMKt;
fINaLp=(1.0./(1.0+KmCaMK./CaMKa));
INaL=c.GNaL.*(V-ENa).*mL.*((1.0-fINaLp).*hL+fINaLp.*hLp);

figure
plot(t,ICaL,'linewidth',2)
hold on
plot(t,IKr,'linewidth',2)
plot(t,IKs,'linewidth',2)
plot(t,INaL,'linewidth',2)
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('time (ms)')
ylabel('current (A/F)')
legend('ICaL','IKr','IKs','INaL')

x1 = find(t==stim_delay); % find interval where the AP begins
x2 = find(floor(V)==floor(V(x1))+3 & t > t(x1)+10,1); % and ends
Area_Ks = trapz(t(x1:x2),IKs(x1:x2));
Area_Kr = trapz(t(x1:x2),IKr(x1:x2));
Area_Ca = trapz(t(x1:x2),ICaL(x1:x2));
Area_NaL = trapz(t(x1:x2),INaL(x1:x2));

IKs_Fraction = Area_Ks/(Area_Ks + Area_Kr);
