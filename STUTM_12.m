% Algoritmo Matlab TG - día
% Estación
% Hora
% Dirección MAC, TI, TF, Duración.

%% 
%Lectura archivo .json
fname = 'rapberry-stutm-D_2021-02-23-export4.json'; % Cambiar 
reportfirebase = jsondecode(fileread(fname)); % Ver al inicio
reporte = reportfirebase; % Ver al final 

% campos de estructura "reporte"
reportefield = fieldnames(reporte);
szreporte = size(reportefield); %Número de estaciones

for i = 1:szreporte(1) % Recorrer cada estación
    %disp(reportefield{i})
    hours = reporte.(reportefield{i});
    hoursfield = fieldnames(hours);
    szhours = size(hoursfield);% Número de horas
    
    tb2 = [];
    for j = 1:szhours(1) % Recorrer cada captura de horas de estación
        %disp (hoursfield{j}) 
        dir_macs = hours.(hoursfield{j});
        %disp(dir_macs);
        macsfield = fieldnames(dir_macs);
        szmacs = size(macsfield) - 1;
               
        for k = 1:szmacs(1) % Recorrer cada dirección MAC
            %disp(macsfield{k})
            MAC = dir_macs.(macsfield{k});
            start_time = MAC.Tiempo_inicial;
            final_time = MAC.Tiempo_final;

            st = datetime(start_time,'Format','HH:mm:ss');
            ft = datetime(final_time,'Format','HH:mm:ss');

            % Llenar srreglos con campos 
            listmacs(k) = string(strrep(macsfield(k),'M_',''));
            listtimes(k,:) = [st,ft];
            listduration(k) = between(st,ft);
            listdurationminutes(k) = minutes(ft - st);
        end

        %nuevos campos en estrucutura hours.hour
        % Hora de captura en string
        newChr = strrep(hoursfield(j),'H_','');  
        hour = strrep(newChr,'_',':');  
        reporte.(reportefield{i}).(hoursfield{j}).horacaptura = string(hour);

        % Tabla con los tiempos de cada mac con duración en cada captura de cada estación 
        Dir_MAC = listmacs';
        Tiempos = listtimes;
        timesort = sort(listtimes);
        Duracion = listduration';
        Duracion_minutos = listdurationminutes';    
        tb1 = table(Dir_MAC, Tiempos, Duracion, Duracion_minutos); 
        tb1 = sortrows(tb1,2); 
        reporte.(reportefield{i}).(hoursfield{j}).registro = tb1;
        
        tb2 = [tb2;tb1];      
        
        % Determinar cantidad de usuarios con los tiempos de entrada y salida
        timeslist = [timesort(:,1); timesort(:,2)];
        timesliststring = string(timeslist);
        timeslist2(1:szmacs) = "EN";
        timeslist2(szmacs+1:szmacs*2) = "SA";
        finaltimelist = [timesliststring timeslist2'];
        finaltimelistsort = sortrows(finaltimelist); 

        numusuarios = [];
        C = 0;
        for z = 1:szmacs*2
            if (finaltimelistsort(z,2) == "EN") 
                num = 1; 
            else 
                num = -1; 
            end 
            numusuarios(z) = C  + num;
            C = numusuarios(z);
        end

        Tiempos_todos = datetime(finaltimelistsort(:,1),'Format','HH:mm:ss');
        Accion = finaltimelistsort(:,2);
        Usuarios = numusuarios';

        tb3 = table(Tiempos_todos, Accion, Usuarios);
        reporte.(reportefield{i}).(hoursfield{j}).registro2  = tb3;  
        
        % Arreglos hora de captura, numero de ususarios y duración promedio 
        % en estación
        hourss(j) = datetime(hour,'Format','HH:mm:ss');
        usershour(j) =  hours.(hoursfield{j}).Usuarios;
        durAve(j) = mean(listdurationminutes); 

        %macs(j) = dirmac;
        
        % Vaciar variables
        clear listmacs listtimes listduration listdurationminutes timeslist2

    end

    %nuevos campos en estrucutura hours
    Hora_captura = hourss';
    Usuarios = usershour';
    Duracion_promedio = durAve';

    tb5 = table(Hora_captura, Usuarios, Duracion_promedio); 
    reporte.(reportefield{i}).registro = tb5; 
    
    reporte.(reportefield{i}).registrocompleto = tb2; 
    
end

reporte.estaciones = string(reportefield);
clearvars -except reporte reportfirebase  

%% Organización registro completo por usuario(MAC)en los tiempos de cada estación

reportefield = fieldnames(reporte);
szreporte = size(reportefield);
szreporte = szreporte(1)-1;

for i = 1:szreporte
    % Revisar si dirección MAC se repite  en el registro completo
    cregister = reporte.(reportefield{i}).registrocompleto;
    [C, ia, ic] = unique(cregister.Dir_MAC,'stable');
    % C = Direcciones MAC sin repetir, ia = indices de primera aparición, 
    % ic = indices de repetición
    
    % Cantidad de veces aparición de cada dirección MAC 
    a_counts = accumarray(ic,1);

    x = 1;
    for j=1:size(ia)
        index = find(ismember(ic,j));% posiciones de aparición
        indexsz = size(index);
        %disp(indexsz(1))
        if(indexsz(1) > 1) % Si hay más de una aparición 
            if (indexsz(1) == 2)% Si hay dos apariciones
                finaltime1 = cregister.Tiempos(index(1),2); 
                starttime2 = cregister.Tiempos(index(2),1);
                difference = minutes(starttime2 - finaltime1);
                %disp(difference)
                if (difference < 5) % Si la diferencia es menor a 5 minutos 
                    %disp(difference)
                    MAC(x) = cregister.Dir_MAC(index(1));
                    ti = cregister.Tiempos(index(1),1);
                    tf = cregister.Tiempos(index(2),2);
                    st(x) = ti;
                    ft(x) = tf;
                    listduration(x) = between(ti,tf);
                    listdurationminutes(x) = minutes(tf - ti);
                    occurrence(x) = 1; % Se considera una misma aparición
                else
                    % Se guardan ambas apariciones
                    %Primero
                    MAC(x) = cregister.Dir_MAC(index(1));
                    ti = cregister.Tiempos(index(1),1);
                    tf = cregister.Tiempos(index(1),2);
                    st(x) = ti;
                    ft(x) = tf;                  
                    listduration(x) = between(ti,tf);
                    listdurationminutes(x) = minutes(tf - ti);
                    occurrence(x) = 1;
                    % Segundo
                    x = x + 1;
                    MAC(x) = cregister.Dir_MAC(index(2));
                    ti = cregister.Tiempos(index(2),1);
                    tf = cregister.Tiempos(index(2),2);
                    st(x) = ti;
                    ft(x) = tf;                    
                    listduration(x) = between(ti,tf);
                    listdurationminutes(x) = minutes(tf - ti);
                    occurrence(x) = 2;
                end
            else % esperar que no pase, mas de 2 casos en un reporte
                %Corregir todos repetidos
                indexsz2 = indexsz(1)-1;
                for z=1:indexsz2
                    finaltime1 = cregister.Tiempos(index(z),2); 
                    starttime2 = cregister.Tiempos(index(z+1),1);
                    difference = minutes(starttime2 - finaltime1);
                    %disp(difference)
                    if (difference < 5) % Si la diferencia es menor a 5 minutos                         
                        if (z == 1)
                            MAC(x) = cregister.Dir_MAC(index(z));
                            ti = cregister.Tiempos(index(z),1);
                            tf = cregister.Tiempos(index(z+1),2);
                            st(x) = ti;
                            ft(x) = tf;                                                     
                            listduration(x) = between(ti,tf);
                            listdurationminutes(x) = minutes(tf - ti);
                            occurrence(x) = 1;                         
                        else
                            x = x - 1;
                            ti = st(x);
                            tf = cregister.Tiempos(index(z+1),2);
                            ft(x) = tf;
                            listduration(x) = between(ti,tf);
                            listdurationminutes(x) = minutes(tf - ti);                            
                            occuren = occurrence(x);
                            occurrence(x) = occuren;                                                               
                        end                    
                    else
                        % Corregir
                        if (z == 1)                            
                            % primero
                            MAC(x) = cregister.Dir_MAC(index(1));
                            ti = cregister.Tiempos(index(z),1);
                            tf = cregister.Tiempos(index(z),2);
                            st(x) = ti;
                            ft(x) = tf;                  
                            listduration(x) = between(ti,tf);
                            listdurationminutes(x) = minutes(tf - ti);
                            occurrence(x) = 1;
                            % Segundo
                            x = x + 1;
                            MAC(x) = cregister.Dir_MAC(index(2));
                            ti = cregister.Tiempos(index(z+1),1);
                            tf = cregister.Tiempos(index(z+1),2);
                            st(x) = ti;
                            ft(x) = tf;                    
                            listduration(x) = between(ti,tf);
                            listdurationminutes(x) = minutes(tf - ti);
                            occurrence(x) = 2;
                        else
                            MAC(x) = cregister.Dir_MAC(index(z+1));
                            ti = cregister.Tiempos(index(z+1),1);
                            tf = cregister.Tiempos(index(z+1),2);
                            st(x) = ti;
                            ft(x) = tf;
                            listduration(x) = between(ti,tf);
                            listdurationminutes(x) = minutes(tf - ti);
                            occuren = occurrence(x-1);
                            occurrence(x) = occuren+1; 
                        end                        
                    end
                    if (z < indexsz2)
                        x = x + 1;
                    end
                end
            end       
        else % Si solo aparece una vez
            MAC(x) = cregister.Dir_MAC(index(1));
            st(x) = cregister.Tiempos(index(1),1);
            ft(x) = cregister.Tiempos(index(1),2);
            %duration(i) = netb.duracion(index(1));
            listduration(x) = cregister.Duracion(index(1));
            listdurationminutes(x) = cregister.Duracion_minutos(index(1));
            occurrence(x) = 1;
        end
        x = x + 1;
    end

    Dir_MAC = MAC'; % Cambiar
    Tiempos = [st' ft'];
    Duracion = listduration';
    Duracion_minutos = listdurationminutes';
    Aparicion = occurrence';
    %Apariciones = a_counts;
    
    % Tabla organizada con registro completo considerando repetición de
    % dirección MAC en reporte 
    tb1 = table(Dir_MAC, Tiempos, Duracion, Duracion_minutos, Aparicion);
    %tb1 = sortrows(tb1,2); 
    reporte.(reportefield{i}).registrocompleto2 = tb1; 
    
    % Tablña organizada con registro de cantidad de usuarios en todos los
    % intervalos con duración promedio
    Duracion_promedio = mean(Duracion_minutos);
    Usuarios = size(Dir_MAC);
    Usuarios = Usuarios(1);
    tb2 = table(Usuarios,Duracion_promedio);
    reporte.(reportefield{i}).registrocompleto3 = tb2; 

end   
    
clearvars -except reporte reportfirebase  

%% Repetición de dirección MAC en las estaciones
% Poder identificar la trazabilidad

reportefield = fieldnames(reporte);
szreporte = size(reportefield); %Número de estaciones

if ((szreporte(1)-1) == 2)
    estacion1 = reporte.Estacion1.registrocompleto2;
    estacion2 = reporte.Estacion2.registrocompleto2;

    [Dir_MAC,ia,ib] = intersect(estacion1.Dir_MAC,estacion2.Dir_MAC,'stable');
    % Dir_MAC = direcciones MAC en ambas estaciones, ia = indices de la
    % aparicion en la primera estacion, ic = indices segunda estación

    Tiempos_Estacion1 = estacion1.Tiempos(ia,:);
    Duracion_Estacion1 = estacion1.Duracion(ia,:); 
    Tiempos_Estacion2 = estacion2.Tiempos(ib,:);
    Duracion_Estacion2 = estacion2.Duracion(ib,:);

    % Tabla con direcciones MAC y tiempos en cada estación 
    tbl = table(Dir_MAC,Tiempos_Estacion1,Tiempos_Estacion2,Duracion_Estacion1,Duracion_Estacion2);
    reporte.registrorepeticion = tbl;

    clearvars -except reporte reportfirebase  
end 
%% Gráficas
%% Graficas opción1
%Mejorar Stem
reportefield = fieldnames(reporte);
szreporte = size(reportefield);
if (szreporte(1) == 2)
    szreporte(1) = 1;
else 
    szreporte(1) = 2;
end

for i = 1:szreporte(1)
    stationname = reporte.estaciones(i);
    stationinfo = reporte.(reportefield{i}).registro;
    
    figure ('Name',stationname,'NumberTitle','off')
    subplot(2,1,1)
    bar(stationinfo.Hora_captura,stationinfo.Usuarios)
    grid on
    ylim([0 max(stationinfo.Usuarios)+1])
    title ("Usuarios captados por intervalo")
    xlabel("Tiempo")
    ylabel("Num usuarios")

    subplot(2,1,2)
    bar(stationinfo.Hora_captura,stationinfo.Duracion_promedio)
    grid on
    ylim([0 max(stationinfo.Duracion_promedio)+1])
    title ("Duración promedio")
    xlabel("Tiempo")
    ylabel("Minutos")

    hours = reporte.(reportefield{i});
    hoursfield = fieldnames(hours);
    szhours = size(hoursfield)-4;% Número de horas
    szhours = szhours(1);
    
    for j = 1:szhours
        hourgrafic = reporte.(reportefield{i}).(hoursfield{j}).horacaptura;

        timeregister = table2array(reporte.(reportefield{i}).(hoursfield{j}).registro2(:,1));
        usersregister = table2array(reporte.(reportefield{i}).(hoursfield{j}).registro2(:,3));

        % figura
        figure('Name',strcat(stationname,{' '},hourgrafic),'NumberTitle','off')
        stairs(timeregister,usersregister)     
        hold on
        stem(timeregister,usersregister,'-og')
        ylim([0 max(usersregister)+1])
        grid on
        title('Registro de usuarios')
        xlabel("Tiempo")
        ylabel("Num usuarios")
    end
end

clearvars -except reporte reportfirebase 
%% Grafica trayectoria, paso de una estación a otra.

reportefield = fieldnames(reporte);
szreporte = size(reportefield); %Número de estaciones

if (szreporte(1) == 4)

    repetition = reporte.registrorepeticion;
    repetitionsz = size(repetition);
    y = [1 1];
    y2 = [2 2];
    for i=1:repetitionsz(1)
        estacion1times = repetition.Tiempos_Estacion1(i,:);
        estacion2times = repetition.Tiempos_Estacion2(i,:);
        %times = [calle45times marlytimes];

        figure1 = figure('Name',repetition.Dir_MAC(i),'NumberTitle','off');
        axes1 = axes('Parent',figure1);
        p1 = line(estacion1times,y);
        p1.Color = 'yellow';    
        hold on 
        p2 = line(estacion2times,y2);
        p2.Color = 'red';
        hold off
        
        ylim([0 3])
        title(repetition.Dir_MAC(i))
        xlabel("Tiempo")
        ylabel("Estaciones")    
        % Set the remaining axes properties
        set(axes1,'YTickLabel',{'','','Est1','','Est2','',''});
        grid on
    end
end
clearvars -except reporte reportfirebase 

%%
% Trayectoria - paso de una estación a otra
% Análisis todas las horas. Hecho, revisar mas de 2 repeticiones de MAC.
% Análisis todas las estaciones. Hecho, solo para 2 estaciones
% Mejorar grafica 1, pendiente
% Grfiac repetición, perndiente


% Documentar todos los codigos
%   AS, Python, Matlab

    