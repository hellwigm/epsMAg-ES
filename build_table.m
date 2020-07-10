clear all, clc

ESvariant{1}        = 'epsMAgES';
 
D = 10;

for jj = 1:1
    filename=ESvariant{jj};
    
    for k=1:28
        strategy    = ESvariant{jj};
        foldername  = ['CEC17_RS_' strategy];
        load([foldername '/' strategy '_RS_on_fun' num2str(k) 'D' num2str(D) '.mat'],'FitT','GB','input')
    
        tab=build_stats(FitT,input);
        
        meanGB = 0;
        for i=1:input.runs
            meanGB=meanGB+GB{i}.evals;
        end
        fevalsUntilGlobalBest(k)=meanGB/input.runs;
        
        Stats(k,:)=[tab(3,:)]; 
    end    

    Statistics = {'Best';'Median';'c';'nu';'Mean';'Std';'Worst';'FR';'vio';'mRTgb'};

    for j=1:28
       nput={
         num2str(Stats(j,1),'%10.5e\n');
        num2str(Stats(j,2),'%10.5e\n');
         ['(' num2str(Stats(j,3)) ',' num2str(Stats(j,4)) ',' num2str(Stats(j,5)) ')'];
         num2str(Stats(j,6),'%10.5e\n');
         num2str(Stats(j,7),'%10.5e\n');
         num2str(Stats(j,8),'%10.5e\n');
         num2str(Stats(j,9),'%10.5e\n');
         num2str(Stats(j,10)*100);
         num2str(Stats(j,11),'%10.5e\n');
         num2str(fevalsUntilGlobalBest(j),'%10.5e\n')}
        if j <= 9
            eval(['C0' num2str(j) '=nput']) 
        else 
            eval(['C' num2str(j) '=nput'])
        end
    end
    
    T=table(C01,C02,C03,C04,C05,C06,C07,C08,C09,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,'Rownames',Statistics)
%    
    filename2=[date '_CEC2017_Results_of_' filename '_in_Dim' num2str(D) '.csv']
    writetable(T,filename2,'WriteRowNames',true)
    
end
    
