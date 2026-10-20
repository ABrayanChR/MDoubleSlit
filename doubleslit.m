function difraccion_optimizada()
    clear; clc; close all;
    
    % --- CONFIGURACIÓN INICIAL ---
    config.MAXIMIZAR_IMAGEN = 1; 
    config.AbsTol = 1e-4; % Tolerancia para integral2
    
    % Parámetros físicos base
    p.lambda = 650e-9;
    p.L = 0.015e-3;
    p.h = 0.02e-3;
    p.d = 0.06e-3;
    p.z = 1;
    p.n = 1;
    p.k = 2*pi*p.n / p.lambda;
    
    % Iniciar Pool Paralelo si no existe (para usar todos los núcleos)
    if isempty(gcp('nocreate'))
        try
            parpool; 
        catch
            disp('No se pudo iniciar parallel pool. Se ejecutará en serie.');
        end
    end

    %% 1. COMPARACIÓN DE MÉTODOS (RS, Fresnel, Fraunhofer)
    disp('--- 1. Calculando Comparación de Métodos ---');
    xpoints = linspace(-0.1, 0.1, 5000);
    
    % A. Rayleigh-Sommerfeld (RS)
    disp('Calculando Rayleigh-Sommerfeld...');
    [I_RS, I0_RS] = calcular_campo_1D(xpoints, p, 'RS', config.AbsTol);
    
    % B. Fresnel
    disp('Calculando Fresnel...');
    [I_Fes, I0_Fes] = calcular_campo_1D(xpoints, p, 'Fresnel', config.AbsTol);
    
    % C. Fraunhofer
    disp('Calculando Fraunhofer...');
    [I_Fra, I0_Fra] = calcular_campo_1D(xpoints, p, 'Fraunhofer', config.AbsTol);

    % --- GRÁFICAS COMPARATIVAS ---
    % 1. RS Absoluta
    titulos = {'Evaluación numérica de la integral de Rayleigh-Sommerfeld', ...
               sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p.z, p.L, p.d, p.lambda)};
    generar_grafica(xpoints, I_RS, [], "$x(m)$", "Intensidad $I$ $(W/m^2)$", titulos, "1evaluacioninicial", config, []);

    % 2. RS Normalizada
    titulos{1} = 'Normalización de la intensidad';
    x_norm_axis = p.d * (xpoints./p.z) / p.lambda;
    generar_grafica(x_norm_axis, I_RS/I0_RS, [], "$\frac{dx}{z \lambda}$", "$I/I_{0}$", titulos, "2intensidadnormalizada", config, []);

    % 3. Comparación Espacial
    titulos{1} = 'Comparación de métodos';
    datos_comp = {I_RS/I0_RS, I_Fes/I0_Fes, I_Fra/I0_Fra};
    legendas = {'Rayleigh-Sommerfeld', 'Aprox. Fresnel', 'Aprox. Fraunhofer'};
    generar_grafica(xpoints, datos_comp, legendas, "$x(m)$", "$I/I_0$", titulos, "3comparacionmetodos_x", config, []);
    
    % 4. Comparación Normalizada
    norm_axis_total = (p.d * xpoints) ./ (p.z * p.lambda);
    generar_grafica(norm_axis_total, datos_comp, legendas, "$\frac{dx}{z \lambda}$", "$I/I_0$", titulos, "3comparacionmetodos_norm", config, []);

    %% 2. VARIANDO Z
    disp('--- 2. Variando Z ---');
    z_factors = [0.001, 1e-4, 1e-6];
    xpoints_high = linspace(-0.2, 0.2, 2000);
    
    for idx = 1:length(z_factors)
        p_local = p;
        p_local.z = p.z * z_factors(idx);
        
        % Ajustar rango X proporcional a Z para ver el patrón
        x_local = xpoints_high * z_factors(idx); 
        
        [I_z, I0_z] = calcular_campo_1D(x_local, p_local, 'RS', config.AbsTol);
        
        titulos = {'Variación de Z', sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p_local.z, p_local.L, p_local.d, p_local.lambda)};
        generar_grafica(x_local, I_z/I0_z, [], "$x(m)$", "$I/I_{0}$", titulos, sprintf("evaluacion_z_%d", idx), config, []);
    end

    %% 3. VARIANDO d
    disp('--- 3. Variando d ---');
    d_factors = [0.001, 1e-4, 1e-6];
    xpoints_d = linspace(-0.1, 0.1, 800);
    
    for idx = 1:length(d_factors)
        p_local = p;
        p_local.d = p.d * d_factors(idx);
        
        % Ajuste de X según el código original
        x_local = xpoints_d * d_factors(idx);
        
        [I_d, I0_d] = calcular_campo_1D(x_local, p_local, 'RS', config.AbsTol);
        
        titulos = {'Variación de d', sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p_local.z, p_local.L, p_local.d, p_local.lambda)};
        generar_grafica(x_local, I_d/I0_d, [], "$x(m)$", "$I/I_{0}$", titulos, sprintf("evaluacion_d_%d", idx), config, []);
    end

    %% 4. VARIANDO L
    disp('--- 4. Variando L ---');
    l_factors = [0.001, 1e-4, 1e-6];
    
    for idx = 1:length(l_factors)
        p_local = p;
        p_local.L = p.L * l_factors(idx);
        
        % Ajuste de X según el código original
        x_local = xpoints * l_factors(idx);
        
        [I_L, I0_L] = calcular_campo_1D(x_local, p_local, 'RS', config.AbsTol);
        
        titulos = {'Variación de L', sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p_local.z, p_local.L, p_local.d, p_local.lambda)};
        generar_grafica(x_local, I_L/I0_L, [], "$x(m)$", "$I/I_{0}$", titulos, sprintf("evaluacion_L_%d", idx), config, []);
    end

    %% 5. VARIANDO LAMBDA
    disp('--- 5. Variando Lambda ---');
    xpoints_lam = linspace(-0.1, 0.1, 200);
    lam_factors = [0.001, 0.01, 0.5];
    
    for idx = 1:length(lam_factors)
        p_local = p;
        p_local.lambda = p.lambda * lam_factors(idx);
        p_local.k = 2*pi*p_local.n / p_local.lambda;
        
        % Nota: Tolerancia relajada a 1e-3 según original para lambda
        [I_lam, I0_lam] = calcular_campo_1D(xpoints_lam, p_local, 'RS', 1e-3);
        
        titulos = {'Variación lambda', sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p_local.z, p_local.L, p_local.d, p_local.lambda)};
        generar_grafica(xpoints_lam, I_lam/I0_lam, [], "$x(m)$", "$I/I_{0}$", titulos, sprintf("evaluacion_lambda_%d", idx), config, []);
    end

    %% 6. REPRESENTACIÓN TRIDIMENSIONAL
    disp('--- 6. Calculando 3D (Optimizado) ---');
    xpoints3D = linspace(-0.05, 0.05, 400); % Resolución original 400x400
    ypoints3D = linspace(-0.05, 0.05, 400);
    
    [Z_int, I0_3D] = calcular_campo_3D(xpoints3D, ypoints3D, p, 'RS', config.AbsTol);
    [X, Y] = meshgrid(xpoints3D, ypoints3D);
    
    figure;
    surf(X, Y, Z_int./I0_3D);
    shading interp;
    titulos = {'Solución de la integral de Rayleigh-Sommerfeld en 3D', ...
               sprintf('(z=%.1e m, L=%.1e m, d=%.1e m, lambda=%.1e m)', p.z, p.L, p.d, p.lambda)};
    title(titulos, 'FontSize', 30);
    xlabel("$x(m)$", 'Interpreter', 'latex', 'FontSize', 26);
    ylabel("$y(m)$", 'Interpreter', 'latex', 'FontSize', 26);
    zlabel("$I/I_{0}$", 'Interpreter', 'latex', 'FontSize', 26);
    colorbar;
    set(gca, 'FontSize', 24);
    if config.MAXIMIZAR_IMAGEN, set(gcf, 'WindowState', 'maximized'); pause(0.5); end
    saveas(gcf, "3intensidadtridimensional", "png");
    
    disp('Proceso completado exitosamente.');
end

%% --- FUNCIONES AUXILIARES Y DE CÁLCULO ---

function [I, I0] = calcular_campo_1D(xpoints, p, metodo, tol)
    % Calcula la intensidad en un array de puntos X (y=0) usando todos los núcleos.
    % Devuelve I (array) e I0 (escalar de referencia en x=0)
    
    % 1. Calcular I0 (en x=0)
    func_kernel0 = obtener_kernel(0, 0, p, metodo);
    val_I0 = integrar_doble_rendija(func_kernel0, p, tol);
    I0 = abs(val_I0)^2;
    
    % 2. Calcular perfil completo (Paralelizado)
    int_vals = zeros(1, length(xpoints));
    y = 0;
    
    % Variables temporales para evitar overhead en parfor
    L = p.L; d = p.d; h = p.h; z = p.z; k = p.k; n = p.n; lambda = p.lambda;
    
    % 'parfor' distribuye las iteraciones entre los núcleos disponibles
    parfor idx = 1:length(xpoints)
        x = xpoints(idx);
        % Reconstruimos el kernel dentro del worker para evitar problemas de serialización
        % Nota: No podemos llamar a 'obtener_kernel' dentro de parfor fácilmente si devuelve func handle complejo
        % así que definimos la lógica aquí o usamos una función pura.
        
        % Lógica directa del Kernel para optimización dentro del loop
        if strcmp(metodo, 'RS')
             r = @(xpr,ypr) sqrt( (x-xpr).^2+(ypr-y).^2+(z).^2 );
             raysom = @(xpr,ypr) (-1/(2.*pi)) .* ( (z.*exp(1i.*k.*n.*r(xpr,ypr)))./r(xpr,ypr).^2 ) .* (1i.*k.*n - 1./r(xpr,ypr));
        elseif strcmp(metodo, 'Fresnel')
             prefactor = (-1i*k*n / (2*pi*z)) * exp(1i*k*n*(z + (x^2 + y^2)/(2*z)));
             kernel_fres = @(xpr, ypr) exp( (1i*k*n/(2*z)) .* (xpr.^2 + ypr.^2 - 2.*(x.*xpr + y.*ypr)) );
             raysom = @(xpr,ypr) prefactor .* kernel_fres(xpr, ypr);
        elseif strcmp(metodo, 'Fraunhofer')
             prefactor = (-1i*k*n / (2*pi*z)) * exp(1i*k*n*(z + (x^2 + y^2)/(2*z)));
             kernel_frau = @(xpr, ypr) exp( ( -1i*k*n/z) .* (x.*xpr + y.*ypr) );
             raysom = @(xpr,ypr) prefactor .* kernel_frau(xpr, ypr);
        else
            raysom = @(xpr,ypr) 0;
        end
        
        % Integración (Doble rendija)
        v1 = integral2(raysom,(-L/2-d/2),(L/2-d/2),(-h/2),(h/2),'AbsTol',tol);
        v2 = integral2(raysom,(-L/2+d/2),(L/2+d/2),(-h/2),(h/2),'AbsTol',tol);
        int_vals(idx) = v1 + v2;
    end
    
    I = abs(int_vals.^2);
end

function [Z_int, I0] = calcular_campo_3D(xvec, yvec, p, metodo, tol)
    % Optimización para matriz 2D
    Nx = length(xvec);
    Ny = length(yvec);
    Z_int = zeros(Ny, Nx);
    
    % Calcular I0
    func_kernel0 = obtener_kernel(0, 0, p, metodo);
    val_I0 = integrar_doble_rendija(func_kernel0, p, tol);
    I0 = abs(val_I0)^2;
    
    L = p.L; d = p.d; h = p.h; z = p.z; k = p.k; n = p.n;
    
    % Usamos parfor en el bucle externo (filas o columnas)
    parfor i = 1:Nx
        x = xvec(i);
        col_vals = zeros(Ny, 1); % Columna temporal
        for j = 1:Ny
            y = yvec(j);
            % Definición local RS (Hardcoded para RS en 3D como en original)
            r = @(xpr,ypr) sqrt( (x-xpr).^2+(ypr-y).^2+(z).^2 );
            raysom = @(xpr,ypr) (-1/(2.*pi)) .* ( (z.*exp(1i.*k.*n.*r(xpr,ypr)))./r(xpr,ypr).^2 ) .* (1i.*k.*n - 1./r(xpr,ypr));
            
            v1 = integral2(raysom,(-L/2-d/2),(L/2-d/2),(-h/2),(h/2),'AbsTol',tol);
            v2 = integral2(raysom,(-L/2+d/2),(L/2+d/2),(-h/2),(h/2),'AbsTol',tol);
            col_vals(j) = abs(v1 + v2)^2;
        end
        Z_int(:, i) = col_vals;
        % fprintf('Columna %d calculada\n', i); % Opcional: ver progreso
    end
end

function func = obtener_kernel(x, y, p, metodo)
    % Genera el function handle para integral2 basado en el método
    k = p.k; n = p.n; z = p.z;
    U = 1; % Amplitud unitaria
    
    if strcmp(metodo, 'RS')
        r = @(xpr,ypr) sqrt( (x-xpr).^2+(ypr-y).^2+(z).^2 );
        func = @(xpr,ypr) (-1/(2.*pi)) .* U .* ...
            ( (z.*exp(1i.*k.*n.*r(xpr,ypr)))./r(xpr,ypr).^2 ) .* ...
            (1i.*k.*n - 1./r(xpr,ypr));
            
    elseif strcmp(metodo, 'Fresnel')
        prefactor = (-1i*k*n / (2*pi*z)) * exp(1i*k*n*(z + (x^2 + y^2)/(2*z)));
        kernel = @(xpr, ypr) exp( (1i*k*n/(2*z)) .* (xpr.^2 + ypr.^2 - 2.*(x.*xpr + y.*ypr)) );
        func = @(xpr,ypr) prefactor .* U .* kernel(xpr, ypr);
        
    elseif strcmp(metodo, 'Fraunhofer')
        prefactor = (-1i*k*n / (2*pi*z)) * exp(1i*k*n*(z + (x^2 + y^2)/(2*z)));
        kernel = @(xpr, ypr) exp( ( -1i*k*n/z) .* (x.*xpr + y.*ypr) );
        func = @(xpr,ypr) prefactor .* U .* kernel(xpr, ypr);
    else
        error('Método no reconocido');
    end
end

function val = integrar_doble_rendija(fun_handle, p, tol)
    % Realiza la integral doble sobre las dos aperturas
    L = p.L; d = p.d; h = p.h;
    % Apertura Izquierda
    val1 = integral2(fun_handle, (-L/2-d/2), (L/2-d/2), (-h/2), (h/2), 'AbsTol', tol);
    % Apertura Derecha
    val2 = integral2(fun_handle, (-L/2+d/2), (L/2+d/2), (-h/2), (h/2), 'AbsTol', tol);
    val = val1 + val2;
end

function generar_grafica(x_data, y_data, leyendas, x_lab, y_lab, titulos, nombre_archivo, config, color_override)
    figure;
    if iscell(y_data) % Si pasamos varias curvas
        hold on
        estilos = {'r', 'b--', 'g:'}; % Estilos para RS, Fresnel, Fraunhofer
        for i = 1:length(y_data)
            plot(x_data, y_data{i}, estilos{i}, 'LineWidth', 2);
        end
        hold off
        if ~isempty(leyendas)
            legend(leyendas, 'Interpreter', 'latex', 'FontSize', 22, 'Location', 'best');
        end
    else
        plot(x_data, y_data, 'LineWidth', 2);
    end
    
    title(titulos, 'Interpreter', 'latex', 'FontSize', 30);
    xlabel(x_lab, 'Interpreter', 'latex', 'FontSize', 26);
    ylabel(y_lab, 'Interpreter', 'latex', 'FontSize', 26);
    grid on; set(gca, 'FontSize', 24);
    
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1) pos(2)*0.8 pos(3) pos(4)*0.85]);
    
    if config.MAXIMIZAR_IMAGEN
        set(gcf, 'WindowState', 'maximized'); 
        pause(0.5); 
    end
    saveas(gcf, nombre_archivo, "png");
end