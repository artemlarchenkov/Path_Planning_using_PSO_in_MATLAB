function sol2 = ParseSolution(sol1, model)
    % PARSESOLUTION Преобразует представление решения и вычисляет свойства пути
    
    % Извлечение координат x и y из решения
    x = sol1.x;
    y = sol1.y;
    
    % Извлечение параметров модели
    xs = model.xs;
    ys = model.ys;
    xt = model.xt;
    yt = model.yt;
    xobs = model.xobs;
    yobs = model.yobs;
    robs = model.robs;
    
    % Формирование массивов координат для построения пути
    XS = [xs, x, xt];
    YS = [ys, y, yt];
    k = numel(XS);
    TS = linspace(0, 1, k);
    
    % Интерполяция для создания более гладкого пути
    tt = linspace(0, 1, 500);
    xx = spline(TS, XS, tt);
    yy = spline(TS, YS, tt);
    
    % Вычисление длины пути
    dx = diff(xx);
    dy = diff(yy);
    L = sum(sqrt(dx.^2 + dy.^2));
    
    % Проверка нарушений ограничений относительно препятствий
    nobs = numel(xobs); % Количество препятствий
    Violation = 0;
    for k = 1:nobs
        d = sqrt((xx - xobs(k)).^2 + (yy - yobs(k)).^2);
        v = max(1 - d / robs(k), 0);
        Violation = Violation + mean(v);
    end
    
    % Формирование структуры свойств пути
    sol2.TS = TS;
    sol2.XS = XS;
    sol2.YS = YS;
    sol2.tt = tt;
    sol2.xx = xx;
    sol2.yy = yy;
    sol2.dx = dx;
    sol2.dy = dy;
    sol2.L = L;
    sol2.Violation = Violation;
    sol2.IsFeasible = (Violation == 0);
end
