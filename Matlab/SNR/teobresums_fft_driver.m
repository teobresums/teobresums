function teobresums_fft_driver
% Generate FFTs
% This scripts generates FFTs from data nd writes files

orange   = [1 0.6 0];
Tau      = [1.0          0.6          0.6             1.0          0.9];
Tau_end  = [0.99            0.99         0.99            0.98         0.96];
alph1    = [0.0035          0.003        0.0025           0.003         0.003];
alph2    = [0.07           0.05         0.04             0.04         0.03];
tapering = [1               1            1               1             1];
npts_end = [800             500          200            200             100];
M        = [2.70            2.60         2.70            2.50           2.50];
clr      = {'b'             'c'          'g'            orange           'r'};
style    = {'-'             '-'          '-'           '-'           '-'};
linewidth= [1                1             1             2              2 ];
ni       = 1;
ne       = 1;
GMsunbyc3= 4.925490947e-6*1e3; % Conversion factor in Hz
Mtot     = 2.7; % Msun
outfile  = 'Data_FT/FT_bns.dat'; 
for n=1
    %fname                        = '../../C++/waveform.dat';
    fname                        = '../../C/bns_q1_s0s0_M2.7_20Hz_postadiab/hlm_interp_l2_m2.txt';
    [o{n}.fM, o{n}.AF, o{n}.Psi] = teobresums_fft(fname,outfile,1,npts_end(n),[Tau(n) Tau_end(n)],[alph1(n) alph2(n)],tapering(n),0);
    o{n}.fMHz                    = o{n}.fM/(M(n).*GMsunbyc3);  
    
    % remove the unphysical region of the FFT
    n0 = 18;
    [~, locs]  = findpeaks(-o{n}.Psi(n0:end));
    
    jmin = locs(1) + n0;
    
    % cut vector
    %%{
    o{n}.fM    =  o{n}.fM(jmin:end);
    o{n}.AF    =  o{n}.AF(jmin:end);
    o{n}.fMHz  =  o{n}.fMHz(jmin:end);
    o{n}.Psi   =  o{n}.Psi(jmin:end);
    o{n}.dPsi  =  Deriv4(o{n}.Psi,o{n}.fM,4);
    o{n}.d2Psi =  Deriv4(o{n}.dPsi,o{n}.fM ,4);
    o{n}.Qomg  =  o{n}.d2Psi.*(o{n}.fM.^2);
    %}
    %------------------------------------------
    % rewriting file adding the kHz column.
    % it overwrites the file written by the FFT
    %------------------------------------------
    %{
    nend = 5500;
    DA_WriteFFTHz(outfile{n},o{n}.fM(1:end-nend),o{n}.Psi(1:end-nend),o{n}.AF(1:end-nend),o{n}.fMHz(1:end-nend))
    %}
end


figure('Name','Phase of the FT in fM')
for n=ni:ne
    semilogx(o{n}.fM,o{n}.Psi,'Color',clr{n},'LineStyle',style{n},'LineWidth',linewidth(n))
hold on;
end
xlabel('$Mf$','FontSize',20,'Interpreter','Latex');
ylabel('$\Psi(fM)$','FontSize',20,'Interpreter','Latex');
xlim([0.0 0.1]);
leg=legend('SLy-135135','SLy-140120','H4-135135','ALF2-140110','MS1b-150100');
set(leg,'Interpreter','Latex');


figure('Name','Amplitude of the FT in fM')
for n=ni:ne
loglog(o{n}.fM,o{n}.AF,'Color',clr{n},'LineStyle',style{n},'LineWidth',linewidth(n))
hold on;
%xlim([0.005 0.1]);
%ylim([0.1 100]);
xlabel('$Mf$','FontSize',20,'Interpreter','Latex');
ylabel('$\tilde{A}(fM)$','FontSize',20,'Interpreter','Latex');
end
leg=legend('SLy-135135','SLy-140120','H4-135135','ALF2-140110','MS1b-150100');
set(leg,'Interpreter','Latex');

end

