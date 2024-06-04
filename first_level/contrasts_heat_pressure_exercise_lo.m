%% Heat Exercise Nacl(Low Intensity Exercise)
%==================

%% Heat Ex Lo SAL
%%-----------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 30
%%--------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat30_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [1 0 0 zeros(1,nConds/2) 0];
    negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 50
%%-------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat50_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [0 1 0 zeros(1,nConds/2) 0];
    negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 70
%% ---------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat70_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [0 0 1 zeros(1,nConds/2) 0];
    negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% Heat Exercise NLX (Low Exercise)
%==================

%% Heat Ex LO NLX
%%------------------

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 30
%%------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat30_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [1 0 0 zeros(1,nConds/2) 0];
    negcondBlock = [-1 0 0 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% Ex low 50
%%%-------------

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat50_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [0 1 0 zeros(1,nConds/2) 0];
    negcondBlock = [0 -1 0 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 70
%% --------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat70_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [0 0 1 zeros(1,nConds/2) 0];
    negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% Heat Exercise NLX > SAL Low Exercise
%%--------------------
% Prep values
tConVec = [];

for r = 1:nSess

    condBlock = [ones(1,nConds/2) zeros(1,nConds/2) 0];
    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negcondBlock = zeroBlock;
    end


    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat_nlx>sal_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Heat Exercise NLX > SAL Low Exercise (VAS 70)
%%--------------------
% Prep values
tConVec = [];

for r = 1:nSess

    condBlock = [0 0 1 zeros(1,nConds/2) 0];
    negcondBlock = [0 0 -1 zeros(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negcondBlock = zeroBlock;
    end


    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'heat70_nlx>sal_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end
%% Parametric Heat SAL (Low exercise)
%==================

% Prep values
tConVec = [];
for r = 1:nSess


    condBlock = [-1 0 1 zeros(1,nConds/2) 0];
    %negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];
    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end



    %name--------------------------------------------------
    tcoName = {'heat_int_sal_lo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Interaction INTENSITY x TREATMENT HEAT (LOW Exercise)
% ===================
% Prep values
tConVec = [];

for r = 1:nSess


    condBlock = [-1 0 1 0 0 0 0];
    negCondBlock = [1 0 -1 0 0 0 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds,1) 0];

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negCondBlock = zeroBlock;
    end


    %name--------------------------------------------------
    tcoName = {'inter_intensity_treatment_heat_lo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end

% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% PRESSURE Ex Lo SAL
%%-----------------
% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Ex low 30

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure30_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 1 0 0 0];
    negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

% Ex low 50

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure50_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 0 1 0 0];
    negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

% Ex low 70

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure70_exLo_sal'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 0 0 1 0];
    negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Pressure Exercise NLX
%==================

% Ex low

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
    negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

% Ex low 30

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure30_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 1 0 0 0];
    negcondBlock = [zeros(1,nConds/2) -1 0 0 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


% Ex low 50

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure50_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 0 1 0 0];
    negcondBlock = [zeros(1,nConds/2) 0 -1 0 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

% Ex low 70

% Prep values
tConVec = [];

subjectConstant = [];%zeros(1,nSess);

params_ses = load(['/projects/crunchie/nold/PEEP/behavioural/MAIN/' sprintf('sub-%02.2d',splitSubs{np}(g)) filesep 'ses-01' filesep sprintf('parameters_sub0%0.02d_bike.mat',splitSubs{np}(g))],'P');

for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure70_exLo_nlx'};

    %convec-------------------------------------------------
    exercise = [params_ses.P.exercise.condition params_ses.P.exercise.condition];
    condBlock = [zeros(1,nConds/2) 0 0 1 0];
    negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Pressure Exercise NLX > SAL Low Exercise
%%--------------------
% Prep values
tConVec = [];

for r = 1:nSess

    condBlock = [zeros(1,nConds/2) ones(1,nConds/2) 0];
    negcondBlock = [zeros(1,nConds/2) -ones(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negcondBlock = zeroBlock;
    end


    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure_nlx>sal_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Pressure Exercise NLX > SAL Low Exercise (VAS 70)
%%--------------------
% Prep values
tConVec = [];

for r = 1:nSess

    condBlock = [zeros(1,nConds/2) 0 0 1 0];
    negcondBlock = [zeros(1,nConds/2) 0 0 -1 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negcondBlock = zeroBlock;
    end


    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    %name--------------------------------------------------
    tcoName = {'pressure70_nlx>sal_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end
%% Parametric Pressure SAL (Low exercise)
%==================

% Prep values
tConVec = [];
for r = 1:nSess


    condBlock = [zeros(1,nConds/2) -1 0 1 0];
    %negcondBlock = [-ones(1,nConds/2) zeros(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds) 0];
    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end



    %name--------------------------------------------------
    tcoName = {'pressure_int_sal_lo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negcondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Interaction INTENSITY x TREATMENT Pressure (LOW Exercise)
% ===================
% Prep values
tConVec = [];

for r = 1:nSess


    condBlock = [0 0 0 -1 0 1 0];
    negCondBlock = [0 0 0 1 0 -1 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds,1) 0];

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negCondBlock = zeroBlock;
    end


    %name--------------------------------------------------
    tcoName = {'inter_intensity_treatment_pressure_lo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end

% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Interaction MODALITY x TREATMENT
%% ===================

% Prep values
tConVec = [];
for r = 1:nSess


    condBlock = [ones(1,nConds/2) -ones(1,nConds/2) 0];
    negCondBlock = [-ones(1,nConds/2) ones(1,nConds/2) 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds,1) 0];

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negCondBlock = zeroBlock;
    end

    %name--------------------------------------------------
    tcoName = {'inter_modality_treatment_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end

% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end

%% Interaction MODALITY x TREATMENT VAS 70
%% ===================

% Prep values
tConVec = [];
for r = 1:nSess


    condBlock = [0 0 1 0 0 -1 0];
    negCondBlock = [0 0 -1 0 0 1 0];
    subjectConstant = [];%zeros(1,nSess);
    zeroBlock = [zeros(1,nConds,1) 0];

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);

    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
        negCondBlock = zeroBlock;
    end

    %name--------------------------------------------------
    tcoName = {'inter_modality_treatment_exLo_VAS70'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec negCondBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end

% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% Heat NLX  > Pressure NLX
%%==================

% Prep values
tConVec = [];


for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);
    condBlock = [ones(1,nConds/2) -ones(1,nConds/2) 0];
    zeroBlock = [zeros(1,nConds) 0];
    subjectConstant = [];%zeros(1,nSess);


    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    %name--------------------------------------------------
    tcoName = {'heat_nlx>pressure_nlx_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end


%% Heat 70 NLX  > Pressure 70 NLX
%%==================

% Prep values
tConVec = [];


for r = 1:nSess

    nMp = n.n_nuis(r);
    mvmnt = zeros(1,nMp);
    condBlock = [0 0 1 0 0 -1 0];
    zeroBlock = [zeros(1,nConds) 0];
    subjectConstant = [];%zeros(1,nSess);


    if exercise(r) == 0
        condBlock;
    elseif exercise(r) == 1
        condBlock = zeroBlock;
    end

    %name--------------------------------------------------
    tcoName = {'heat70_nlx>pressure70_nlx_exLo'};

    %convec------------------------------------------------

    if treatment_order(splitSubs{np}(g)) == 1
        if r < floor(nSess/2)+1
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        end

    elseif treatment_order(splitSubs{np}(g)) == 0
        if r < floor(nSess/2)+1
            tConVec = [tConVec zeroBlock mvmnt subjectConstant];
        elseif r > floor(nSess/2)
            tConVec = [tConVec condBlock mvmnt subjectConstant];
        end
    end

end

%fill t-cons into template
for co = 1:size(tcoName,2)
    tco = tco + 1;
    template.spm.stats.con.consess{tco}.tcon.name    = tcoName{:};
    template.spm.stats.con.consess{tco}.tcon.convec  = tConVec;
    template.spm.stats.con.consess{tco}.tcon.sessrep = 'none';
end


% Sanity Check: Contrast Vector length
if length(tConVec) ~= length(condBlock)*8 + sum(n.n_nuis(1:end))
    error('Contrast Length is not coherent');
end