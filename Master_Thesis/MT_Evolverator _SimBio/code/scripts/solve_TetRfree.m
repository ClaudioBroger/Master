syms lexAtot lexAon est Kdest
sol = solve((lexAtot - lexAon)*(est - lexAon) == Kdest*lexAon,lexAon);
lexA = logspace(0,3,10);
Kdest = 0.9;
lexAtot = logspace(0,3,10);
solLexA = subs(sol(1),'lexAtot',lexAtot);
solLexA = subs(solLexA,'Kdest',Kdest);
solLexA1 = subs(solLexA,'est',1);
solLexA10 = subs(solLexA,'est',10);
solLexA100 = subs(solLexA,'est',100);
solLexA1000 = subs(solLexA,'est',1000);

close all
plot(log10(lexAtot),log10(solLexA1));
hold on; plot(log10(lexAtot),log10(solLexA10))
hold on; plot(log10(lexAtot),log10(solLexA100))
hold on; plot(log10(lexAtot),log10(solLexA1000))
xlabel('log10([LexA])')
ylabel('log10([LexAon])')
legend('est = 1','est = 10','est = 100','est = 1000')

%%
syms aTc TetR TetRon Kd
solve(TetRon*(aTc - TetR + TetRon) == Kd*(TetR - TetRon),TetRon)
sol = solve(TetRon*(aTc - TetR + TetRon) == Kd*(TetR - TetRon),TetRon);
TetR = logspace(0,3,10);
aTc = 10;
Kd = 15;
solTetR = subs(sol(2),'TetR',TetR);
solTetR = subs(solTetR,'Kd',Kd);
solTetR1 = subs(solTetR,'aTc',1);
solTetR10 = subs(solTetR,'aTc',10);
solTetR100 = subs(solTetR,'aTc',100);
solTetR1000 = subs(solTetR,'aTc',1000);
plot(log10(TetR),log10(solTetR1));
hold on; plot(log10(TetR),log10(solTetR10))
hold on; plot(log10(TetR),log10(solTetR100));
hold on; plot(log10(TetR),log10(solTetR1000));
legend('aTc = 1','aTc = 10','aTc = 100','aTc = 1000')
xlabel('log10([TetR])')
ylabel('log10([TetRon])')