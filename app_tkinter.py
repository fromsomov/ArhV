"""bobrinsky_analyzer_v6_scientific.py
Анализатор сосудов Bobrinsky с научной классификацией Цетлина
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, Menu, colorchooser
import threading
import queue
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('TkAgg')
from datetime import datetime
import pathlib
import re
import logging
from collections import Counter

# Попробуем импортировать tkinterdnd2 для drag-and-drop
try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    HAVE_DND = True
except ImportError:
    HAVE_DND = False
    print("Библиотека tkinterdnd2 не установлена. Drag-and-drop не будет работать.")
    print("Установите её: pip install tkinterdnd2")

# Попробуем импортировать numpy-stl для экспорта 3D моделей
try:
    from stl import mesh
    HAVE_STL = True
except ImportError:
    HAVE_STL = False
    print("Библиотека numpy-stl не установлена. Экспорт в STL не будет работать.")
    print("Установите её: pip install numpy-stl")

# Импорты для DXF обработки
import ezdxf
from scipy.interpolate import interp1d, CubicSpline
from scipy.integrate import simpson, trapezoid
import math
from collections import defaultdict

# ============================================================================
# КЛАССИФИКАЦИЯ ЦЕТЛИНА (научная шкала)
# ============================================================================

TSETLIN_CLASSIFICATION_L = [
    {
        'group': 'I',
        'start_l': 0.024,
        'center_l': 0.035,
        'end_l': 0.049,
        'quality_name': 'Супермалый-2',
        'mobility_class': '1 – «супермалые» (менее 0,097 л)',
        'description': 'Сосуды для хранения ароматических веществ'
    },
    {
        'group': 'II',
        'start_l': 0.049,
        'center_l': 0.071,
        'end_l': 0.097,
        'quality_name': 'Супермалый-1',
        'mobility_class': '1 – «супермалые» (менее 0,097 л)',
        'description': 'Сосуды для хранения ароматических веществ'
    },
    {
        'group': 'III',
        'start_l': 0.097,
        'center_l': 0.137,
        'end_l': 0.194,
        'quality_name': 'Очень очень малый',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'IV',
        'start_l': 0.194,
        'center_l': 0.274,
        'end_l': 0.389,
        'quality_name': 'Очень малый',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'V',
        'start_l': 0.389,
        'center_l': 0.552,
        'end_l': 0.782,
        'quality_name': 'Малый',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'VI',
        'start_l': 0.782,
        'center_l': 1.105,
        'end_l': 1.565,
        'quality_name': 'Мало-средний',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'VII',
        'start_l': 1.565,
        'center_l': 2.210,
        'end_l': 3.125,
        'quality_name': 'Средний-1',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'VIII',
        'start_l': 3.125,
        'center_l': 4.420,
        'end_l': 6.250,
        'quality_name': 'Средний-2',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'IX',
        'start_l': 6.250,
        'center_l': 8.840,
        'end_l': 12.500,
        'quality_name': 'Средний-3',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'X',
        'start_l': 12.500,
        'center_l': 17.680,
        'end_l': 25.000,
        'quality_name': 'Средний-4',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'XI',
        'start_l': 25.000,
        'center_l': 35.360,
        'end_l': 50.0,
        'quality_name': 'Больше-средний',
        'mobility_class': '2 – «мобильные» (0,097 л – 50,0 л)',
        'description': 'Легко перемещаются одним взрослым человеком'
    },
    {
        'group': 'XII',
        'start_l': 50.0,
        'center_l': 70.7,
        'end_l': 100.0,
        'quality_name': 'Большой',
        'mobility_class': '3 – «ограниченно-мобильные» (50,0 л – 200,0 л)',
        'description': 'Требуют усилий минимум двух человек'
    },
    {
        'group': 'XIII',
        'start_l': 100.0,
        'center_l': 141.4,
        'end_l': 200.0,
        'quality_name': 'Очень большой',
        'mobility_class': '3 – «ограниченно-мобильные» (50,0 л – 200,0 л)',
        'description': 'Требуют усилий минимум двух человек'
    },
    {
        'group': 'XIV',
        'start_l': 200.0,
        'center_l': 282.9,
        'end_l': 400.0,
        'quality_name': 'Очень очень большой',
        'mobility_class': '4 – «мало-мобильные» (200,0 л – 800,0 л)',
        'description': 'Перемещались крайне редко, только пустыми'
    },
    {
        'group': 'XV',
        'start_l': 400.0,
        'center_l': 565.8,
        'end_l': 800.0,
        'quality_name': 'Гигантский',
        'mobility_class': '4 – «мало-мобильные» (200,0 л – 800,0 л)',
        'description': 'Перемещались крайне редко, только пустыми'
    },
    {
        'group': 'XVI',
        'start_l': 800.0,
        'center_l': 1131.5,
        'end_l': 1600.0,
        'quality_name': 'Супер-1',
        'mobility_class': '5 – «условно-мобильные» (800,0 л – 3200,0 л)',
        'description': 'Перемещаются только в незаполненном виде'
    },
    {
        'group': 'XVII',
        'start_l': 1600.0,
        'center_l': 2263.0,
        'end_l': 3200.0,
        'quality_name': 'Супер-2',
        'mobility_class': '5 – «условно-мобильные» (800,0 л – 3200,0 л)',
        'description': 'Перемещаются только в незаполненном виде'
    },
    {
        'group': 'XVIII',
        'start_l': 3200.0,
        'center_l': 4526.0,
        'end_l': 6400.0,
        'quality_name': 'Сверх-1',
        'mobility_class': '6 – «стационарные» (3200,0 л – 25000,0 л)',
        'description': 'В принципе не предполагают перемещения'
    },
    {
        'group': 'XIX',
        'start_l': 6400.0,
        'center_l': 9052.0,
        'end_l': 12800.0,
        'quality_name': 'Сверх-2',
        'mobility_class': '6 – «стационарные» (3200,0 л – 25000,0 л)',
        'description': 'В принципе не предполагают перемещения'
    },
    {
        'group': 'XX',
        'start_l': 12800.0,
        'center_l': 18104.0,
        'end_l': 25000.0,
        'quality_name': 'Сверх-3',
        'mobility_class': '6 – «стационарные» (3200,0 л – 25000,0 л)',
        'description': 'В принципе не предполагают перемещения'
    }
]

# ============================================================================
# ЦВЕТОВАЯ СХЕМА В СОВРЕМЕННОМ СТИЛЕ
# ============================================================================

MODERN_PALETTE = {
    'primary': '#2c3e50',
    'primary_light': '#3498db',
    'primary_dark': '#1a252f',
    'secondary': '#95a5a6',
    'accent': '#e74c3c',
    'success': '#27ae60',
    'warning': '#f39c12',
    'danger': '#c0392b',
    'light': '#ecf0f1',
    'dark': '#2c3e50',
    'bg_light': '#ffffff',
    'bg_card': '#f8f9fa',
    'border': '#dee2e6'
}

GRADIENT = ['#3498db', '#2980b9', '#1f639b', '#154a7d', '#0c355f']

# ============================================================================
# КЛАСС ПРОФИЛЬНОЙ ГРУППЫ
# ============================================================================

class ProfileGroup:
    """Класс для группировки профилей сосудов"""
    
    def __init__(self, name):
        self.name = name
        self.profiles = []
        self.files = []
        
    def add_profile(self, profile, file_path):
        self.profiles.append(profile)
        self.files.append(file_path)
        
    def remove_profile(self, file_path):
        if file_path in self.files:
            idx = self.files.index(file_path)
            self.files.pop(idx)
            return self.profiles.pop(idx)
        return None
        
    def get_stats(self):
        if not self.profiles:
            return {}
        
        volumes = [p.get('volume', 0) for p in self.profiles if p]
        heights = [np.max(p.get('y', [0])) for p in self.profiles if p]
        
        return {
            'count': len(self.profiles),
            'avg_volume': np.mean(volumes) if volumes else 0,
            'avg_height': np.mean(heights) if heights else 0,
            'total_volume': sum(volumes) if volumes else 0
        }

# ============================================================================
# КОРРЕКТНЫЙ РАСЧЁТ ОБЪЁМА (ИСПРАВЛЕННЫЙ)
# ============================================================================

class CorrectVolumeCalculator:
    """Исправленный калькулятор объёма с разными методами"""
    
    def __init__(self, y_coords, r_coords):
        self.y = np.asarray(y_coords, dtype=np.float64)
        self.r = np.asarray(r_coords, dtype=np.float64)
        
        # Убираем нулевые и отрицательные радиусы
        self.r = np.maximum(self.r, 0.001)
        
        # Сортируем по высоте
        sort_idx = np.argsort(self.y)
        self.y = self.y[sort_idx]
        self.r = self.r[sort_idx]
        
        # Убираем дубликаты высот
        unique_y, unique_idx = np.unique(self.y, return_index=True)
        self.y = unique_y
        self.r = self.r[unique_idx]
        
        # Проверяем данные
        if len(self.y) < 2:
            raise ValueError("Слишком мало точек для расчета объема")
        
        # Убедимся, что высота начинается с 0
        if self.y[0] > 0.01:
            self.y = np.insert(self.y, 0, 0.0)
            self.r = np.insert(self.r, 0, self.r[0])
        
        # Инициализируем интерполяторы
        self._init_interpolators()
        
        # Диагностика
        print(f"\nИнициализация калькулятора:")
        print(f"  Количество точек: {len(self.y)}")
        print(f"  Высота: от {self.y[0]:.2f} до {self.y[-1]:.2f} см")
        print(f"  Радиус: от {self.r[0]:.2f} до {self.r[-1]:.2f} см")
    
    def _init_interpolators(self):
        """Инициализация интерполяторов с проверкой"""
        try:
            # Для сплайна нужно минимум 4 точки
            if len(self.y) >= 4:
                self.spline = CubicSpline(self.y, self.r, bc_type='natural')
                print("  Используется кубический сплайн")
            else:
                self.spline = interp1d(self.y, self.r, kind='cubic', 
                                     fill_value='extrapolate', bounds_error=False)
                print("  Используется кубическая интерполяция (мало точек)")
        except Exception as e:
            print(f"  Ошибка создания сплайна: {e}")
            print("  Используется линейная интерполяция")
            self.spline = interp1d(self.y, self.r, kind='linear', 
                                 fill_value='extrapolate', bounds_error=False)
        
        # Линейный интерполятор для некоторых методов
        self.linear = interp1d(self.y, self.r, kind='linear', 
                             fill_value='extrapolate', bounds_error=False)
    
    def method_disks(self, y_max=None):
        """
        Метод дисков - интегрирование площадей кругов
        Использует исходные точки без интерполяции
        """
        if y_max is None:
            y_max = self.y[-1]
        
        # Находим индекс последней точки до y_max
        mask = self.y <= y_max
        y_slice = self.y[mask]
        r_slice = self.r[mask]
        
        if len(y_slice) < 2:
            return 0.0
        
        # Если y_max не совпадает с последней точкой, добавляем её
        if y_max > y_slice[-1] and y_max < self.y[-1]:
            # Интерполируем радиус для y_max
            r_max = float(self.linear(y_max))
            y_slice = np.append(y_slice, y_max)
            r_slice = np.append(r_slice, r_max)
        
        # Площади дисков
        areas = np.pi * r_slice**2
        
        # Интегрирование методом трапеций по исходным точкам
        return trapezoid(areas, y_slice)
    
    def method_frustums(self, y_max=None):
        """
        Метод усеченных конусов
        Объем между двумя сечениями: V = π/3 * h * (r1² + r1*r2 + r2²)
        """
        if y_max is None:
            y_max = self.y[-1]
        
        # Находим индекс последней точки до y_max
        mask = self.y <= y_max
        y_slice = self.y[mask]
        r_slice = self.r[mask]
        
        if len(y_slice) < 2:
            return 0.0
        
        # Если y_max не совпадает с последней точкой
        if y_max > y_slice[-1]:
            r_max = float(self.linear(y_max))
            y_slice = np.append(y_slice, y_max)
            r_slice = np.append(r_slice, r_max)
        
        volume = 0.0
        for i in range(len(y_slice) - 1):
            h = y_slice[i + 1] - y_slice[i]
            if h <= 0:
                continue
            
            r1 = r_slice[i]
            r2 = r_slice[i + 1]
            volume += (np.pi / 3.0) * h * (r1**2 + r1 * r2 + r2**2)
        
        return volume
    
    def method_trapezoidal(self, y_max=None, n_points=2000):
        """
        Метод трапеций с интерполяцией
        Большое количество точек для точности
        """
        if y_max is None:
            y_max = self.y[-1]
        
        # Создаем равномерную сетку с большим количеством точек
        y_fine = np.linspace(0, y_max, n_points)
        r_fine = self.spline(y_fine)
        r_fine = np.maximum(r_fine, 0.0)  # Защита от отрицательных радиусов
        
        # Площади сечений
        areas = np.pi * r_fine**2
        
        # Интегрирование методом трапеций
        return trapezoid(areas, y_fine)
    
    def method_simpson(self, y_max=None, n_points=501):
        """
        Метод Симпсона (парабол)
        Нечетное количество точек для точности
        """
        if y_max is None:
            y_max = self.y[-1]
        
        # Для метода Симпсона желательно нечетное число точек
        if n_points % 2 == 0:
            n_points += 1
        
        y_fine = np.linspace(0, y_max, n_points)
        r_fine = self.spline(y_fine)
        r_fine = np.maximum(r_fine, 0.0)
        
        # Площади сечений
        areas = np.pi * r_fine**2
        
        # Интегрирование методом Симпсона
        return simpson(areas, y_fine)
    
    def method_spline_integral(self, y_max=None, n_points=1001):
        """
        Интеграл сплайна с аналитическим интегрированием
        Самый точный метод для гладких профилей
        """
        if y_max is None:
            y_max = self.y[-1]
        
        # Для точности используем много точек
        y_fine = np.linspace(0, y_max, n_points)
        r_fine = self.spline(y_fine)
        r_fine = np.maximum(r_fine, 0.0)
        
        # Площади сечений
        areas = np.pi * r_fine**2
        
        # Используем метод Симпсона для максимальной точности
        return simpson(areas, y_fine)
    
    def calculate_all_methods(self, y_max=None):
        """
        Вычисление объема всеми методами с диагностикой
        """
        if y_max is None:
            y_max = self.y[-1]
        
        print(f"\nРАСЧЕТ ОБЪЕМА ДО ВЫСОТЫ {y_max:.2f} см:")
        print("-" * 60)
        
        methods = {
            'disks': ('Метод дисков (исходные точки)', self.method_disks),
            'frustums': ('Метод усеченных конусов', self.method_frustums),
            'trapezoidal': ('Метод трапеций (2000 точек)', 
                          lambda y: self.method_trapezoidal(y, 2000)),
            'simpson': ('Метод Симпсона (501 точка)', 
                       lambda y: self.method_simpson(y, 501)),
            'spline': ('Интеграл сплайна (1001 точка, рекоменд.)', 
                      lambda y: self.method_spline_integral(y, 1001))
        }
        
        results = {}
        
        # Сначала вычисляем все объемы
        for name, (description, method) in methods.items():
            try:
                results[name] = method(y_max)
            except Exception as e:
                results[name] = None
                print(f"❌ Ошибка в методе '{description}': {e}")
        
        # Выводим результаты
        for name, (description, _) in methods.items():
            if results[name] is not None:
                print(f"✅ {description}:")
                print(f"   Объем: {results[name]/1000:.6f} л")
                print(f"   Объем: {results[name]:.2f} см³")
        
        # Сравниваем с эталонным методом (сплайн)
        if results['spline'] is not None:
            print(f"\n📊 СРАВНЕНИЕ С МЕТОДОМ СПЛАЙНА:")
            print("-" * 40)
            
            for name in ['disks', 'frustums', 'trapezoidal', 'simpson']:
                if results[name] is not None:
                    diff = results[name] - results['spline']
                    diff_percent = (diff / results['spline']) * 100
                    diff_abs = abs(diff_percent)
                    
                    if diff_abs < 0.1:
                        status = "✓ Очень близко"
                    elif diff_abs < 1.0:
                        status = "✓ Близко"
                    elif diff_abs < 5.0:
                        status = "⚠ Приемлемо"
                    else:
                        status = "⚠ Заметная разница"
                    
                    print(f"{methods[name][0]}:")
                    print(f"  Разница: {diff:+.2f} см³ ({diff_percent:+.3f}%)")
                    print(f"  Статус: {status}")
        
        return results
    
    def calculate_volume(self, method_name, y_max=None):
        """Вычисление объема указанным методом с FIX для методов"""
        if y_max is None:
            y_max = self.y[-1]
        
        # ВАЖНОЕ ИСПРАВЛЕНИЕ: Используем правильные методы
        if method_name == 'disks':
            return self.method_disks(y_max)
        elif method_name == 'frustums':
            return self.method_frustums(y_max)
        elif method_name == 'trapezoidal':
            return self.method_trapezoidal(y_max, n_points=2000)
        elif method_name == 'simpson':
            return self.method_simpson(y_max, n_points=501)
        elif method_name == 'spline':
            return self.method_spline_integral(y_max, n_points=1001)
        else:
            # По умолчанию используем метод дисков
            return self.method_disks(y_max)

# ============================================================================
# КЛАСС ДЛЯ КАСТОМНОЙ ПАНЕЛИ ИНСТРУМЕНТОВ
# ============================================================================

class CustomNavigationToolbar:
    """Простая панель инструментов только с кнопками сохранения"""
    
    def __init__(self, canvas, parent):
        self.canvas = canvas
        self.toolbar = ttk.Frame(parent)
        
        save_btn = ttk.Button(self.toolbar, text="💾 Сохранить", 
                             command=self.save_figure)
        save_btn.pack(side=tk.LEFT, padx=2)
        
        self.toolbar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def save_figure(self):
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("All files", "*.*")
            ],
            title="Сохранить график"
        )
        
        if filename:
            try:
                self.canvas.figure.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Успех", f"График сохранён в {filename}")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить: {str(e)}")

# ============================================================================
# ГЛАВНЫЙ КЛАСС ПРИЛОЖЕНИЯ
# ============================================================================

class BobrinskyAnalyzer:
    """Анализатор сосудов Bobrinsky с научной классификацией Цетлина"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Bobrinsky - Анализатор сосудов v6.0 (с классификацией Цетлина)")
        self.root.geometry("1600x900")
        
        # Настройка современной цветовой схемы
        self.setup_modern_style()
        
        # Инициализация данных
        self.profiles = {}
        self.groups = {}  # Убрана группа "Без группы" по умолчанию
        self.current_profile = None
        self.current_group = None
        self.volume_calculator = None
        
        # Состояние раскрытия групп
        self.expanded_groups = set()
        
        # 3D данные для экспорта
        self.X_surface = None
        self.Y_surface = None
        self.Z_surface = None
        
        # Настройки производительности и 3D
        self.settings = {
            'rdp_epsilon': 0.02,  # Параметр упрощения RDP
            'min_profile_points': 50,  # Минимальное количество точек профиля
            'max_profile_points': 500,  # Максимальное количество точек (для 3D)
            '3d_segments': 30,  # Количество сегментов в 3D-модели
            'enable_3d_optimization': True,  # Включить оптимизацию 3D
        }
        
        # Параметры 3D-визуализации
        self.alpha_3d_var = tk.DoubleVar(value=0.8)
        self.surface_color_hex = '#3498db'  # Цвет по умолчанию
        self.surface_style_3d_var = tk.StringVar(value='solid')  # 'solid' или 'wireframe'
        self.projection_type_3d_var = tk.StringVar(value='persp')  # 'persp' или 'ortho'
        self.segments_y_var = tk.IntVar(value=30)
        self.segments_theta_var = tk.IntVar(value=30)
        self.density_var = tk.IntVar(value=2)
        self.show_axes_3d_var = tk.BooleanVar(value=True)  # Включение/выключение осей 3D
        
        # Классификация Цетлина
        self.tsetlin_classification = TSETLIN_CLASSIFICATION_L
        
        # Очередь для обработки
        self.processing_queue = queue.Queue()
        
        # КРИТИЧЕСКАЯ ОШИБКА: В исходном коде метод по умолчанию был 'spline',
        # но он переопределялся в разных местах. Фиксируем это:
        self.method_var = tk.StringVar(value="spline")
        
        # Создание интерфейса
        self.create_interface()
        
        # Запуск обработки очереди
        self.start_queue_processor()
        
        # Настройка drag-and-drop если доступно
        if HAVE_DND:
            self.setup_drag_drop()
    
    def roman_to_int(self, roman):
        """Преобразование римских чисел в арабские"""
        roman_numerals = {
            'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5,
            'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10,
            'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15,
            'XVI': 16, 'XVII': 17, 'XVIII': 18, 'XIX': 19, 'XX': 20
        }
        return roman_numerals.get(roman, 0)
    
    def get_tsetlin_classification(self, volume_cm3):
        """Определение качественной группы объема по классификации Цетлина"""
        volume_l = volume_cm3 / 1000.0
        
        for class_data in self.tsetlin_classification:
            if class_data['start_l'] <= volume_l <= class_data['end_l']:
                # Проверяем, находится ли объем в диапазоне строгого качества
                # (ближе к центру, чем к краям интервала)
                center = class_data['center_l']
                start = class_data['start_l']
                end = class_data['end_l']
                interval_length = end - start
                distance_to_center = abs(volume_l - center)
                
                # Если расстояние до центра меньше 25% интервала, считаем строгим качеством
                is_strict_quality = distance_to_center < (interval_length * 0.25)
                
                return {
                    'group': class_data['group'],
                    'group_name': class_data['quality_name'],
                    'volume_l': volume_l,
                    'center_l': center,
                    'start_l': start,
                    'end_l': end,
                    'is_strict_quality': is_strict_quality,
                    'mobility_class': class_data['mobility_class'],
                    'description': class_data['description']
                }
        
        # Если объем вне диапазонов (очень маленький или очень большой)
        if volume_l < self.tsetlin_classification[0]['start_l']:
            class_data = self.tsetlin_classification[0]
            return {
                'group': class_data['group'],
                'group_name': f"{class_data['quality_name']} (ниже диапазона)",
                'volume_l': volume_l,
                'center_l': class_data['center_l'],
                'start_l': class_data['start_l'],
                'end_l': class_data['end_l'],
                'is_strict_quality': False,
                'mobility_class': class_data['mobility_class'],
                'description': class_data['description']
            }
        else:  # volume_l > self.tsetlin_classification[-1]['end_l']
            class_data = self.tsetlin_classification[-1]
            return {
                'group': class_data['group'],
                'group_name': f"{class_data['quality_name']} (выше диапазона)",
                'volume_l': volume_l,
                'center_l': class_data['center_l'],
                'start_l': class_data['start_l'],
                'end_l': class_data['end_l'],
                'is_strict_quality': False,
                'mobility_class': class_data['mobility_class'],
                'description': class_data['description']
            }
    
    def set_axes_equal(self, ax):
        """Исправленная функция: Устанавливает равный масштаб для осей x, y, z"""
        try:
            # Получаем текущие лимиты
            if hasattr(self, 'X_surface') and self.X_surface is not None:
                # Центрируем сосуд по X и Z
                x_min, x_max = np.min(self.X_surface), np.max(self.X_surface)
                z_min, z_max = np.min(self.Z_surface), np.max(self.Z_surface)
                y_min, y_max = np.min(self.Y_surface), np.max(self.Y_surface)
                
                # Для симметричных сосудов X и Z должны быть центрированы в 0
                # Берем максимальное абсолютное значение для X и Z
                max_xz = max(abs(x_min), abs(x_max), abs(z_min), abs(z_max))
                x_limits = (-max_xz, max_xz)
                z_limits = (-max_xz, max_xz)
                y_limits = (y_min, y_max)
            else:
                x_limits = ax.get_xlim3d()
                y_limits = ax.get_ylim3d()
                z_limits = ax.get_zlim3d()
            
            x_range = abs(x_limits[1] - x_limits[0])
            y_range = abs(y_limits[1] - y_limits[0])
            z_range = abs(z_limits[1] - z_limits[0])
            
            # Находим максимальный диапазон среди ВСЕХ трех осей
            max_range = max(x_range, y_range, z_range)
            
            # Вычисляем средние точки
            # X и Z центрированы в 0, Y центрируется по своей середине
            x_middle = 0.0
            z_middle = 0.0
            y_middle = np.mean(y_limits)
            
            # Устанавливаем новые лимиты с ОДИНАКОВЫМ диапазоном для всех осей
            ax.set_xlim3d([x_middle - max_range/2, x_middle + max_range/2])
            ax.set_ylim3d([z_middle - max_range/2, z_middle + max_range/2])  # y -> z
            ax.set_zlim3d([y_middle - max_range/2, y_middle + max_range/2])  # z -> y
            
            # Устанавливаем равное соотношение сторон
            ax.set_box_aspect([1, 1, 1])
            
        except Exception as e:
            logging.error(f"Ошибка при установке равных осей: {e}")
            try:
                # Резервный метод
                if hasattr(self, 'X_surface') and self.X_surface is not None:
                    # Центрируем сосуд по X и Z
                    x_min, x_max = np.min(self.X_surface), np.max(self.X_surface)
                    z_min, z_max = np.min(self.Z_surface), np.max(self.Z_surface)
                    y_min, y_max = np.min(self.Y_surface), np.max(self.Y_surface)
                    
                    # Для симметричных сосудов X и Z должны быть центрированы в 0
                    max_xz = max(abs(x_min), abs(x_max), abs(z_min), abs(z_max))
                    
                    max_range = max(max_xz * 2, y_max - y_min)
                    
                    y_middle = (y_min + y_max) / 2
                    
                    ax.set_xlim3d([-max_range/2, max_range/2])
                    ax.set_ylim3d([-max_range/2, max_range/2])
                    ax.set_zlim3d([y_middle - max_range/2, y_middle + max_range/2])
                    ax.set_box_aspect([1, 1, 1])
            except:
                pass
    
    def setup_drag_drop(self):
        """Настройка drag-and-drop из проводника (только если tkinterdnd2 доступен)"""
        try:
            # Регистрируем главное окно как цель для перетаскивания файлов
            self.root.drop_target_register(DND_FILES)
            self.root.dnd_bind('<<Drop>>', self.handle_drop)
            
            # Также регистрируем дерево
            self.tree.drop_target_register(DND_FILES)
            self.tree.dnd_bind('<<Drop>>', self.handle_tree_drop)
            
            self.status_var.set("Drag-and-drop включен. Перетащите DXF файлы!")
        except Exception as e:
            print(f"Ошибка настройки drag-and-drop: {e}")
            self.status_var.set("Drag-and-drop не доступен")
    
    def handle_drop(self, event):
        """Обработка перетаскивания файлов в главное окно"""
        try:
            # Получаем список файлов из события
            files = self.parse_dropped_files(event.data)
            self.add_files_to_current_group(files)
        except Exception as e:
            print(f"Ошибка обработки перетаскивания: {e}")
            messagebox.showerror("Ошибка", f"Не удалось обработать файлы: {str(e)}")
    
    def handle_tree_drop(self, event):
        """Обработка перетаскивания файлов в дерево"""
        try:
            item = self.tree.identify_row(event.y)
            files = self.parse_dropped_files(event.data)
            
            if item:
                item_data = self.tree.item(item)
                tags = item_data.get('tags', [])
                
                if 'group' in tags:
                    # Перетаскиваем в группу
                    group_name = item_data['text']
                    self.add_files_to_group(files, group_name)
                elif 'file' in tags:
                    # Перетаскиваем в файл - добавляем в группу файла
                    group_id = self.tree.parent(item)
                    if group_id:
                        group_name = self.tree.item(group_id, 'text')
                        self.add_files_to_group(files, group_name)
            else:
                # Перетаскиваем в пустое место
                self.add_files_to_current_group(files)
        except Exception as e:
            print(f"Ошибка обработки перетаскивания в дерево: {e}")
    
    def parse_dropped_files(self, data):
        """Парсинг строки с файлами из события drag-and-drop"""
        # tkinterdnd2 возвращает строку с файлами, разделенными пробелами
        # Файлы с пробелами в путях заключены в фигурные скобки
        
        files = []
        i = 0
        n = len(data)
        
        while i < n:
            if data[i] == '{':
                # Находим закрывающую скобку
                j = data.find('}', i + 1)
                if j == -1:
                    # Нет закрывающей скобки, берем остаток
                    files.append(data[i:].strip())
                    break
                files.append(data[i+1:j].strip())
                i = j + 1
            else:
                # Находим следующий пробел
                j = data.find(' ', i)
                if j == -1:
                    # Последний файл
                    files.append(data[i:].strip())
                    break
                files.append(data[i:j].strip())
                i = j
            
            # Пропускаем пробелы
            while i < n and data[i] == ' ':
                i += 1
        
        return [os.path.normpath(f) for f in files if f]
    
    def add_files_to_group(self, files, group_name):
        """Добавить файлы в указанную группу"""
        added_count = 0
        for file_path in files:
            if file_path.lower().endswith('.dxf'):
                file_path = os.path.normpath(file_path)
                if file_path not in self.profiles:
                    if group_name not in self.groups:
                        self.groups[group_name] = ProfileGroup(group_name)
                    self.groups[group_name].add_profile(None, file_path)
                    self.profiles[file_path] = None
                    added_count += 1
        
        self.update_tree()
        if added_count > 0:
            self.status_var.set(f"Добавлено {added_count} файлов в группу '{group_name}'")
    
    def add_files_to_current_group(self, files):
        """Добавить файлы в текущую или новую группу"""
        if not files:
            return
        
        if self.current_group and self.current_group in self.groups:
            self.add_files_to_group(files, self.current_group)
        else:
            # Создаем новую группу с именем из первой папки
            first_file = files[0]
            folder_name = os.path.basename(os.path.dirname(first_file))
            if not folder_name or folder_name == '.':
                folder_name = "Новая группа"
            
            base_name = folder_name
            counter = 1
            while base_name in self.groups:
                base_name = f"{folder_name}_{counter}"
                counter += 1
            
            self.current_group = base_name
            self.groups[base_name] = ProfileGroup(base_name)
            self.add_files_to_group(files, base_name)
    
    def simplify_profile_rdp(self, points, epsilon=0.01):
        """
        Упрощение профиля с помощью алгоритма Рамера-Дугласа-Пьюкера.
        Уменьшает количество точек без существенной потери точности.
        
        points: массив точек (N, 2) с координатами [x, y]
        epsilon: допуск упрощения (чем больше, тем сильнее упрощение)
        """
        if len(points) < 3:
            return points
        
        # Находим точку с максимальным расстоянием
        dmax = 0
        index = 0
        start, end = points[0], points[-1]
        
        for i in range(1, len(points)-1):
            d = self.perpendicular_distance(points[i], start, end)
            if d > dmax:
                index = i
                dmax = d
        
        # Рекурсивно упрощаем
        if dmax > epsilon:
            left = self.simplify_profile_rdp(points[:index+1], epsilon)
            right = self.simplify_profile_rdp(points[index:], epsilon)
            return np.vstack((left[:-1], right))
        else:
            return np.array([start, end])
    
    def perpendicular_distance(self, point, line_start, line_end):
        """Вычисление перпендикулярного расстояния от точки до линии"""
        x, y = point
        x1, y1 = line_start
        x2, y2 = line_end
        
        if x1 == x2 and y1 == y2:
            return np.sqrt((x - x1)**2 + (y - y1)**2)
        
        # Формула расстояния от точки до линии
        return np.abs((x2-x1)*(y1-y) - (x1-x)*(y2-y1)) / np.sqrt((x2-x1)**2 + (y2-y1)**2)
    
    def setup_modern_style(self):
        style = ttk.Style()
        style.configure('Modern.TFrame', background=MODERN_PALETTE['bg_light'])
        style.configure('Modern.TLabel', background=MODERN_PALETTE['bg_light'], 
                       foreground=MODERN_PALETTE['dark'], font=('Segoe UI', 10))
        style.configure('Modern.TButton', font=('Segoe UI', 10), padding=8)
        style.configure('Modern.TNotebook', background='white')
        style.configure('Modern.TNotebook.Tab', 
                       padding=[20, 8],
                       font=('Segoe UI', 11, 'bold'),
                       background=MODERN_PALETTE['light'],
                       foreground=MODERN_PALETTE['dark'])
        style.map('Modern.TNotebook.Tab',
                 background=[('selected', 'white')],
                 foreground=[('selected', MODERN_PALETTE['primary_dark'])])
        style.configure('Card.TFrame', background=MODERN_PALETTE['bg_card'], 
                       relief='flat', borderwidth=1)
        style.configure('Card.TLabel', background=MODERN_PALETTE['bg_card'],
                       font=('Segoe UI', 11, 'bold'), foreground=MODERN_PALETTE['primary_dark'])
        style.configure('Treeview', font=('Segoe UI', 10), rowheight=28)
        style.configure('Treeview.Heading', font=('Segoe UI', 11, 'bold'),
                       background=MODERN_PALETTE['primary_light'], 
                       foreground='white')
        
        plt.rcParams.update({
            'axes.prop_cycle': plt.cycler('color', GRADIENT),
            'axes.facecolor': '#FFFFFF',
            'figure.facecolor': '#FFFFFF',
            'axes.edgecolor': MODERN_PALETTE['primary'],
            'axes.labelcolor': MODERN_PALETTE['primary_dark'],
            'text.color': MODERN_PALETTE['dark'],
            'xtick.color': MODERN_PALETTE['primary_dark'],
            'ytick.color': MODERN_PALETTE['primary_dark'],
            'grid.color': '#D6DBDF',
            'grid.alpha': 0.5
        })
    
    def create_interface(self):
        # Главный контейнер
        main_container = ttk.Frame(self.root, style='Modern.TFrame')
        main_container.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Верхняя панель инструментов
        self.create_toolbar(main_container)
        
        # Рабочая область
        workspace = ttk.Frame(main_container)
        workspace.pack(fill=tk.BOTH, expand=True, pady=(10, 0))
        
        # Левая панель - группы и профили
        left_panel = self.create_left_panel(workspace)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # Центральная область с вкладками (РЕОРГАНИЗОВАНО)
        center_panel = self.create_center_panel(workspace)
        center_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Статус бар
        self.create_status_bar(main_container)
    
    def create_toolbar(self, parent):
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill=tk.X, pady=(0, 10))
        
        buttons = [
            ("📁 Добавить DXF", self.add_dxf_files, MODERN_PALETTE['primary']),
            ("⚙️ Обработать", self.process_files, MODERN_PALETTE['warning']),
            ("📊 Экспорт Excel", self.export_excel, MODERN_PALETTE['success']),
            ("🔄 Обновить", self.update_plots, MODERN_PALETTE['secondary']),
            ("⚡ Производительность", self.show_performance_settings, MODERN_PALETTE['accent']),
            ("🔬 Цетлин", self.show_tsetlin_info, '#8e44ad'),
            ("❓ Справка", self.show_help, MODERN_PALETTE['primary_light']),
            ("🧪 Тест объемов", self.test_volume_calculation, '#9b59b6')  # Новая кнопка
        ]
        
        for text, command, color in buttons:
            btn = tk.Button(toolbar, text=text, command=command,
                          bg=color, fg='white',
                          font=('Segoe UI', 10, 'bold'),
                          relief='flat', padx=15, pady=8,
                          cursor='hand2', bd=0,
                          activebackground=self.lighten_color(color, 0.2),
                          activeforeground='white')
            btn.pack(side=tk.LEFT, padx=3)
    
    def create_left_panel(self, parent):
        panel = ttk.Frame(parent, width=320)
        
        # Заголовок
        title_frame = ttk.Frame(panel)
        title_frame.pack(fill=tk.X, pady=(0, 10))
        
        title = ttk.Label(title_frame, text="📂 Группы профилей", 
                         font=('Segoe UI', 12, 'bold'),
                         foreground=MODERN_PALETTE['primary'])
        title.pack(side=tk.LEFT)
        
        # Дерево групп и профилей
        tree_frame = ttk.Frame(panel)
        tree_frame.pack(fill=tk.BOTH, expand=True)
        
        self.tree = ttk.Treeview(tree_frame, show='tree', selectmode='extended')
        
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Привязка событий
        self.tree.bind('<<TreeviewSelect>>', self.on_tree_select)
        self.tree.bind('<Double-1>', self.on_tree_double_click)
        
        # Контекстное меню для дерева
        self.tree_menu = Menu(self.root, tearoff=0)
        self.tree_menu.add_command(label="Создать новую группу", command=self.create_new_group)
        self.tree_menu.add_command(label="Добавить файлы в группу", command=self.add_dxf_files)
        self.tree_menu.add_separator()
        self.tree_menu.add_command(label="Сортировать по имени", command=self.sort_groups_by_name)
        self.tree_menu.add_command(label="Переместить в другую группу", command=self.move_to_group)
        self.tree_menu.add_command(label="Удалить выбранное", command=self.delete_selected)
        self.tree_menu.add_separator()
        self.tree_menu.add_command(label="Обновить дерево", command=self.update_tree)
        self.tree.bind('<Button-3>', self.show_tree_menu)
        
        # Контекстное меню для групп
        self.group_menu = Menu(self.root, tearoff=0)
        self.group_menu.add_command(label="Удалить группу", command=self.delete_group)
        self.group_menu.add_command(label="Переименовать группу", command=self.rename_group)
        self.group_menu.add_separator()
        self.group_menu.add_command(label="Добавить файлы", command=self.add_dxf_files)
        
        return panel
    
    def create_center_panel(self, parent):
        """Создание центральной панели с РЕОРГАНИЗОВАННОЙ структурой вкладок"""
        panel = ttk.Frame(parent)
        
        self.notebook = ttk.Notebook(panel, style='Modern.TNotebook')
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Вкладка 1: ОБЪЁМ (переименована с "Профиль")
        self.tab_volume = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_volume, text="📊 Объём")
        self.setup_volume_tab()  # Новая структура с подвкладками
        
        # Вкладка 2: 3D Модель (без изменений)
        self.tab_3d = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_3d, text="🏺 3D Модель")
        self.setup_3d_tab()
        
        # Вкладка 3: Морфология (без изменений)
        self.tab_morphology = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_morphology, text="🔬 Морфология")
        self.setup_morphology_tab()
        
        return panel
    
    def setup_volume_tab(self):
        """Настройка вкладки 'Объём' с подвкладками"""
        # Создаем Notebook для подвкладок внутри вкладки "Объём"
        volume_notebook = ttk.Notebook(self.tab_volume)
        volume_notebook.pack(fill=tk.BOTH, expand=True)
        
        # Подвкладка 1: Профиль (График и Объемы)
        self.subtab_profile = ttk.Frame(volume_notebook)
        volume_notebook.add(self.subtab_profile, text="📐 Профиль и Объемы")
        self.setup_profile_subtab()
        
        # Подвкладка 2: Таблицы (перенесено из старых "Результатов")
        self.subtab_tables = ttk.Frame(volume_notebook)
        volume_notebook.add(self.subtab_tables, text="📋 Таблицы")
        self.setup_tables_subtab()
        
        # Подвкладка 3: Графики (перенесено из старых "Результатов", с изменениями)
        self.subtab_charts = ttk.Frame(volume_notebook)
        volume_notebook.add(self.subtab_charts, text="📈 Графики")
        self.setup_charts_subtab()
        
        # Подвкладка 4: Классификация Цетлина (НОВАЯ)
        self.subtab_tsetlin = ttk.Frame(volume_notebook)
        volume_notebook.add(self.subtab_tsetlin, text="🎯 Классификация Цетлина")
        self.setup_tsetlin_subtab()
    
    def setup_profile_subtab(self):
        """Настройка подвкладки 'Профиль и Объемы' (старая вкладка 'Профиль')"""
        # Создаем PanedWindow для разделения с поддержкой изменения размеров
        paned = ttk.PanedWindow(self.subtab_profile, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Левая панель - график профиля
        left_frame = ttk.Frame(paned)
        paned.add(left_frame, weight=3)
        
        # Создаём фигуру matplotlib для профиля
        self.fig_profile = Figure(figsize=(10, 6), dpi=100)
        self.ax_profile = self.fig_profile.add_subplot(111)
        self.ax_profile.set_facecolor('#FFFFFF')
        self.ax_profile.grid(True, alpha=0.3, color='#D6DBDF')
        
        self.canvas_profile = FigureCanvasTkAgg(self.fig_profile, left_frame)
        self.canvas_profile.draw()
        
        # Кастомная панель инструментов
        self.toolbar_profile = CustomNavigationToolbar(self.canvas_profile, left_frame)
        self.canvas_profile.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Правая панель - управление объёмом
        right_frame = ttk.Frame(paned)
        paned.add(right_frame, weight=1)
        
        self.setup_volume_panel_in_profile(right_frame)
        
        # Привязка событий для перетаскивания линии уровня
        self.canvas_profile.mpl_connect('button_press_event', self.on_profile_click)
        self.canvas_profile.mpl_connect('motion_notify_event', self.on_profile_drag)
        self.canvas_profile.mpl_connect('button_release_event', self.on_profile_release)
        
        self.dragging_level = False
    
    def setup_volume_panel_in_profile(self, parent):
        """Панель управления объемом в подвкладке 'Профиль'"""
        # Создаём контейнер с прокруткой
        canvas = tk.Canvas(parent, highlightthickness=0, bg='white')
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Заполняем прокручиваемую область
        self.fill_volume_controls(scrollable_frame)
    
    def setup_tables_subtab(self):
        """Настройка подвкладки 'Таблицы' (перенесено из 'Результатов')"""
        # Создаем Notebook для таблиц внутри подвкладки
        tables_notebook = ttk.Notebook(self.subtab_tables)
        tables_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Таблица 1: Результаты анализа
        tab_results = ttk.Frame(tables_notebook)
        tables_notebook.add(tab_results, text="📊 Результаты анализа")
        self.setup_results_table(tab_results)
        
        # Таблица 2: Классификация Цетлина
        tab_tsetlin_table = ttk.Frame(tables_notebook)
        tables_notebook.add(tab_tsetlin_table, text="🎯 Шкала Цетлина")
        self.setup_tsetlin_table(tab_tsetlin_table)
    
    def setup_results_table(self, parent):
        """Таблица результатов анализа (из старой вкладки 'Результаты')"""
        columns = ['Профиль', 'Группа', 'Объём (л)', 'Объём (см³)', 'Высота', 'Диаметр', 'Метод', 'Группа Цетлина']
        
        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.results_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', height=20)
        
        col_widths = [180, 100, 80, 90, 70, 70, 100, 120]
        for col, width in zip(columns, col_widths):
            self.results_tree.heading(col, text=col)
            self.results_tree.column(col, width=width)
        
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.results_tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.results_tree.xview)
        self.results_tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        hsb.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.tree_menu_results = Menu(self.root, tearoff=0)
        self.tree_menu_results.add_command(label="Копировать", command=self.copy_tree_selection)
        self.results_tree.bind('<Button-3>', self.show_tree_menu_results)
    
    def setup_tsetlin_table(self, parent):
        """Таблица с полной шкалой классификации Цетлина"""
        columns = ['Группа', 'Начало (л)', 'Центр (л)', 'Конец (л)', 'Качество', 'Класс мобильности']
        
        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.tsetlin_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', height=20)
        
        col_widths = [70, 100, 100, 100, 120, 200]
        for col, width in zip(columns, col_widths):
            self.tsetlin_tree.heading(col, text=col)
            self.tsetlin_tree.column(col, width=width)
        
        # Заполняем таблицу данными классификации
        for class_data in self.tsetlin_classification:
            self.tsetlin_tree.insert('', 'end', values=(
                class_data['group'],
                f"{class_data['start_l']:.3f}",
                f"{class_data['center_l']:.3f}",
                f"{class_data['end_l']:.3f}",
                class_data['quality_name'],
                class_data['mobility_class']
            ))
        
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tsetlin_tree.yview)
        self.tsetlin_tree.configure(yscrollcommand=vsb.set)
        
        self.tsetlin_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Добавляем информацию о текущем сосуде, если есть
        if self.current_profile and 'tsetlin_classification' in self.current_profile:
            tsetlin_info = self.current_profile['tsetlin_classification']
            info_frame = ttk.Frame(parent)
            info_frame.pack(fill=tk.X, pady=(10, 0), padx=10)
            
            info_text = f"Текущий сосуд: Группа {tsetlin_info['group']} ({tsetlin_info['group_name']}), Объём: {tsetlin_info['volume_l']:.3f} л"
            ttk.Label(info_frame, text=info_text, font=('Segoe UI', 10, 'bold'),
                     foreground=MODERN_PALETTE['primary']).pack()
    
    def setup_charts_subtab(self):
        """Настройка подвкладки 'Графики' (с улучшениями)"""
        self.fig_charts = Figure(figsize=(12, 8), dpi=100)
        
        # 4 графика: улучшенные версии
        self.ax_chart1 = self.fig_charts.add_subplot(221)
        self.ax_chart2 = self.fig_charts.add_subplot(222)
        self.ax_chart3 = self.fig_charts.add_subplot(223)
        self.ax_chart4 = self.fig_charts.add_subplot(224)
        
        self.canvas_charts = FigureCanvasTkAgg(self.fig_charts, self.subtab_charts)
        self.canvas_charts.draw()
        
        self.toolbar_charts = CustomNavigationToolbar(self.canvas_charts, self.subtab_charts)
        self.canvas_charts.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        btn_frame = ttk.Frame(self.subtab_charts)
        btn_frame.pack(fill=tk.X, pady=(0, 10))
        ttk.Button(btn_frame, text="🔄 Обновить графики", 
                  command=self.update_results_charts).pack()
    
    def setup_tsetlin_subtab(self):
        """Настройка подвкладки 'Классификация Цетлина'"""
        # Создаем контейнер с прокруткой
        canvas = tk.Canvas(self.subtab_tsetlin, highlightthickness=0, bg='white')
        scrollbar = ttk.Scrollbar(self.subtab_tsetlin, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Заполняем прокручиваемую область
        self.fill_tsetlin_info(scrollable_frame)
    
    def fill_tsetlin_info(self, parent):
        """Заполнение информации о классификации Цетлина"""
        # Заголовок
        title_frame = ttk.Frame(parent)
        title_frame.pack(fill=tk.X, pady=(10, 20), padx=10)
        
        ttk.Label(title_frame, text="🎯 Научная классификация сосудов по Ю.Б. Цетлину", 
                 font=('Segoe UI', 16, 'bold'),
                 foreground=MODERN_PALETTE['primary']).pack()
        
        # Текущая классификация (если есть сосуд)
        if self.current_profile and 'tsetlin_classification' in self.current_profile:
            tsetlin_info = self.current_profile['tsetlin_classification']
            
            current_card = self.create_card(parent, "📊 Текущая качественная группа")
            
            # Отображаем информацию о классификации
            info_text = f"""
            Группа качества: {tsetlin_info['group']} ({tsetlin_info['group_name']})
            Объём сосуда: {tsetlin_info['volume_l']:.3f} л ({tsetlin_info['volume_l']*1000:.1f} см³)
            Диапазон группы: {tsetlin_info['start_l']:.3f} – {tsetlin_info['end_l']:.3f} л
            Центр качества: {tsetlin_info['center_l']:.3f} л
            
            Класс мобильности: {tsetlin_info['mobility_class']}
            """
            
            if tsetlin_info['is_strict_quality']:
                info_text += "\n✅ Объём соответствует строгому качеству (близок к центру группы)"
            else:
                info_text += "\n⚠️ Объём находится в переходной зоне между группами"
            
            info_text += f"\n\n📝 {tsetlin_info['description']}"
            
            ttk.Label(current_card, text=info_text, justify=tk.LEFT,
                     font=('Segoe UI', 11)).pack(pady=10)
        
        # Полная шкала классификации
        scale_card = self.create_card(parent, "📏 Полная шкала качественных групп")
        
        # Создаем таблицу с группами
        for class_data in self.tsetlin_classification:
            group_frame = ttk.Frame(scale_card)
            group_frame.pack(fill=tk.X, pady=3, padx=5)
            
            # Цветная метка группы
            color_label = tk.Label(group_frame, text="  ", bg=self.get_tsetlin_color(class_data['group']))
            color_label.pack(side=tk.LEFT, padx=(0, 10))
            
            # Информация о группе
            group_info = f"Группа {class_data['group']}: {class_data['quality_name']} " \
                        f"({class_data['start_l']:.3f} – {class_data['end_l']:.3f} л)"
            ttk.Label(group_frame, text=group_info, font=('Segoe UI', 10)).pack(side=tk.LEFT)
        
        # Классы мобильности
        mobility_card = self.create_card(parent, "🚶 Классы мобильности сосудов")
        
        mobility_classes = [
            ("1 – «супермалые» (менее 0,097 л)", "Сосуды для хранения ароматических веществ"),
            ("2 – «мобильные» (0,097 л – 50,0 л)", "Легко перемещаются одним взрослым человеком в заполненном состоянии"),
            ("3 – «ограниченно-мобильные» (50,0 л – 200,0 л)", "Для перемещения требуются усилия минимум двух человек"),
            ("4 – «мало-мобильные» (200,0 л – 800,0 л)", "Перемещались крайне редко, только пустыми"),
            ("5 – «условно-мобильные» (800,0 л – 3200,0 л)", "Перемещаются только в незаполненном виде усилиями нескольких человек"),
            ("6 – «стационарные» (3200,0 л – 25000,0 л)", "В принципе не предполагают перемещения")
        ]
        
        for i, (cls, desc) in enumerate(mobility_classes):
            cls_frame = ttk.Frame(mobility_card)
            cls_frame.pack(fill=tk.X, pady=2, padx=5)
            
            ttk.Label(cls_frame, text=cls, font=('Segoe UI', 10, 'bold'),
                     foreground=MODERN_PALETTE['primary']).pack(anchor='w')
            ttk.Label(cls_frame, text=desc, font=('Segoe UI', 9)).pack(anchor='w', padx=10)
    
    def get_tsetlin_color(self, group):
        """Получение цвета для группы Цетлина"""
        colors = ['#3498db', '#2ecc71', '#e74c3c', '#f39c12', '#9b59b6', '#1abc9c',
                 '#34495e', '#e67e22', '#16a085', '#8e44ad', '#27ae60', '#d35400',
                 '#c0392b', '#2980b9', '#f1c40f', '#7f8c8d', '#95a5a6', '#2c3e50']
        
        try:
            group_num = int(group) - 1
            return colors[group_num % len(colors)]
        except:
            return '#95a5a6'
    
    def setup_3d_tab(self):
        """Вкладка 3D модели с элементами управления - ИСПРАВЛЕННАЯ ВЕРСИЯ"""
        # Создаем основной контейнер с разделением
        main_paned = ttk.PanedWindow(self.tab_3d, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # ЛЕВАЯ панель - график (теперь график слева)
        plot_panel = ttk.Frame(main_paned)
        main_paned.add(plot_panel, weight=4)  # Больший вес для графика
        
        # ПРАВАЯ панель - элементы управления (теперь справа)
        control_panel = ttk.Frame(main_paned, width=280)
        main_paned.add(control_panel, weight=1)
        
        # Настройка графика ПЕРВОЙ (левая часть)
        self.setup_3d_plot_area(plot_panel)
        
        # Настройка панели управления ВТОРОЙ (правая часть)
        self.setup_3d_control_panel(control_panel)
    
    def setup_3d_control_panel(self, parent):
        """Создание панели управления для 3D-визуализации с ВКЛЮЧЕНИЕМ/ВЫКЛЮЧЕНИЕМ ОСЕЙ"""
        # Создаём контейнер с прокруткой
        canvas = tk.Canvas(parent, highlightthickness=0, bg='white')
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Заголовок панели управления
        title_label = ttk.Label(scrollable_frame, text="⚙️ Управление 3D-видом",
                               font=('Segoe UI', 12, 'bold'),
                               foreground=MODERN_PALETTE['primary'])
        title_label.pack(pady=(10, 15))
        
        # Карточка: Визуализация
        viz_card = self.create_card(scrollable_frame, "🎨 Визуализация")
        
        # Прозрачность
        alpha_frame = ttk.Frame(viz_card)
        alpha_frame.pack(fill=tk.X, pady=5)
        ttk.Label(alpha_frame, text="Прозрачность:").pack(side=tk.LEFT)
        self.alpha_3d_var = tk.DoubleVar(value=0.8)
        alpha_slider = ttk.Scale(alpha_frame, from_=0.1, to=1.0, 
                               orient=tk.HORIZONTAL, variable=self.alpha_3d_var,
                               length=150, command=lambda v: self.update_3d_plot())
        alpha_slider.pack(side=tk.RIGHT, padx=5)
        alpha_value = ttk.Label(alpha_frame, text="0.8")
        alpha_value.pack(side=tk.RIGHT)
        
        # Обновление метки при изменении слайдера
        def update_alpha_label(v):
            alpha_value.config(text=f"{float(v):.1f}")
        self.alpha_3d_var.trace('w', lambda *args: update_alpha_label(self.alpha_3d_var.get()))
        
        # Стиль отображения
        style_frame = ttk.Frame(viz_card)
        style_frame.pack(fill=tk.X, pady=5)
        ttk.Label(style_frame, text="Стиль:").pack(side=tk.LEFT)
        self.surface_style_3d_var = tk.StringVar(value='solid')
        style_combo = ttk.Combobox(style_frame, textvariable=self.surface_style_3d_var,
                                  values=['solid', 'wireframe'], state='readonly', width=15)
        style_combo.pack(side=tk.RIGHT, padx=5)
        style_combo.bind('<<ComboboxSelected>>', lambda e: self.update_3d_plot())
        
        # Выбор цвета
        color_frame = ttk.Frame(viz_card)
        color_frame.pack(fill=tk.X, pady=5)
        ttk.Label(color_frame, text="Цвет:").pack(side=tk.LEFT)
        self.color_button = tk.Button(color_frame, text="Выбрать", 
                                     command=self.choose_3d_color,
                                     bg=self.surface_color_hex, fg='white',
                                     font=('Segoe UI', 9))
        self.color_button.pack(side=tk.RIGHT, padx=5)
        
        # Включение/выключение осей (НОВОЕ)
        axes_frame = ttk.Frame(viz_card)
        axes_frame.pack(fill=tk.X, pady=5)
        ttk.Checkbutton(axes_frame, text="Показывать оси координат",
                       variable=self.show_axes_3d_var,
                       command=self.update_3d_plot).pack(side=tk.LEFT)
        
        # Карточка: Сетка и сегментация
        grid_card = self.create_card(scrollable_frame, "📐 Сетка и сегментация")
        
        # Сегменты по высоте
        segments_frame = ttk.Frame(grid_card)
        segments_frame.pack(fill=tk.X, pady=5)
        ttk.Label(segments_frame, text="Сегментов по высоте:").pack(side=tk.LEFT)
        self.segments_y_var = tk.IntVar(value=30)
        segments_spin = ttk.Spinbox(segments_frame, from_=10, to=200,
                                  textvariable=self.segments_y_var, width=10)
        segments_spin.pack(side=tk.RIGHT, padx=5)
        segments_spin.bind('<Return>', lambda e: self.update_3d_plot())
        segments_spin.bind('<FocusOut>', lambda e: self.update_3d_plot())
        
        # Сегменты по окружности
        theta_frame = ttk.Frame(grid_card)
        theta_frame.pack(fill=tk.X, pady=5)
        ttk.Label(theta_frame, text="Сегментов по окружности:").pack(side=tk.LEFT)
        self.segments_theta_var = tk.IntVar(value=30)
        theta_spin = ttk.Spinbox(theta_frame, from_=10, to=100,
                               textvariable=self.segments_theta_var, width=10)
        theta_spin.pack(side=tk.RIGHT, padx=5)
        theta_spin.bind('<Return>', lambda e: self.update_3d_plot())
        theta_spin.bind('<FocusOut>', lambda e: self.update_3d_plot())
        
        # Плотность сетки
        density_frame = ttk.Frame(grid_card)
        density_frame.pack(fill=tk.X, pady=5)
        ttk.Label(density_frame, text="Плотность сетки:").pack(side=tk.LEFT)
        self.density_var = tk.IntVar(value=2)
        density_slider = ttk.Scale(density_frame, from_=1, to=10,
                                 orient=tk.HORIZONTAL, variable=self.density_var,
                                 length=150, command=lambda v: self.update_3d_plot())
        density_slider.pack(side=tk.RIGHT, padx=5)
        
        # Карточка: Камера и проекция
        camera_card = self.create_card(scrollable_frame, "📷 Камера и проекция")
        
        # Проекция
        projection_frame = ttk.Frame(camera_card)
        projection_frame.pack(fill=tk.X, pady=5)
        ttk.Label(projection_frame, text="Проекция:").pack(side=tk.LEFT)
        self.projection_type_3d_var = tk.StringVar(value='persp')
        projection_combo = ttk.Combobox(projection_frame, textvariable=self.projection_type_3d_var,
                                       values=['persp', 'ortho'], state='readonly', width=15)
        projection_combo.pack(side=tk.RIGHT, padx=5)
        projection_combo.bind('<<ComboboxSelected>>', lambda e: self.update_3d_plot())
        
        # Кнопки управления камерой
        camera_buttons_frame = ttk.Frame(camera_card)
        camera_buttons_frame.pack(fill=tk.X, pady=10)
        
        ttk.Button(camera_buttons_frame, text="🔄 Сброс вида",
                  command=self.reset_3d_view).pack(fill=tk.X, pady=2)
        ttk.Button(camera_buttons_frame, text="📐 Изометрический вид",
                  command=self.set_isometric_view).pack(fill=tk.X, pady=2)
        ttk.Button(camera_buttons_frame, text="👆 Вид сверху",
                  command=self.set_top_view).pack(fill=tk.X, pady=2)
        
        # Карточка: Экспорт
        export_card = self.create_card(scrollable_frame, "💾 Экспорт")
        
        ttk.Button(export_card, text="📸 Сохранить снимок",
                  command=self.save_3d_snapshot).pack(fill=tk.X, pady=3)
        
        # Кнопка экспорта STL
        if HAVE_STL:
            ttk.Button(export_card, text="📦 Экспорт STL",
                      command=self.export_3d_model).pack(fill=tk.X, pady=3)
        else:
            btn = ttk.Button(export_card, text="📦 Экспорт STL (установите numpy-stl)",
                           command=lambda: messagebox.showerror("Ошибка", 
                           "Для экспорта в STL установите библиотеку: pip install numpy-stl"))
            btn.pack(fill=tk.X, pady=3)
        
        # Информация о модели (упрощенная версия)
        info_card = self.create_card(scrollable_frame, "ℹ️ Информация")
        
        self.model_info_label = ttk.Label(info_card, text="Модель не загружена",
                                         wraplength=250, justify=tk.LEFT)
        self.model_info_label.pack(fill=tk.X, pady=5)
    
    def setup_3d_plot_area(self, parent):
        """Настройка области для 3D графика - МАКСИМАЛЬНО РАСШИРЕННАЯ"""
        # Создаём фигуру matplotlib для 3D
        self.fig_3d = Figure(figsize=(10, 8), dpi=100)
        self.ax_3d = self.fig_3d.add_subplot(111, projection='3d')
        
        # Настройка 3D графика
        self.ax_3d.set_facecolor('#FFFFFF')
        
        # Создаем canvas и toolbar
        self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, parent)
        self.toolbar_3d = CustomNavigationToolbar(self.canvas_3d, parent)
        
        # Упаковываем canvas чтобы он занимал все доступное пространство
        self.canvas_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
    
    def setup_morphology_tab(self):
        """Вкладка морфологии"""
        morphology_frame = ttk.Frame(self.tab_morphology)
        morphology_frame.pack(fill=tk.BOTH, expand=True, padx=20, pady=20)
        
        title = ttk.Label(morphology_frame, text="🔬 Анализ морфологии сосудов", 
                         font=('Segoe UI', 16, 'bold'),
                         foreground=MODERN_PALETTE['primary'])
        title.pack(pady=(0, 20))
        
        info_text = """Функции анализа морфологии:

• Анализ кривизны профиля
• Определение типа сосуда
• Сравнение с эталонными формами
• Статистический анализ геометрии
• Классификация по морфологическим признакам

Данный модуль находится в стадии разработки
и будет доступен в следующем обновлении."""
        
        info_label = ttk.Label(morphology_frame, text=info_text,
                              font=('Segoe UI', 11),
                              justify=tk.LEFT,
                              background='white',
                              padding=20)
        info_label.pack(fill=tk.BOTH, expand=True)
    
    def fill_volume_controls(self, parent):
        # Карточка: Метод расчёта
        method_card = self.create_card(parent, "📏 Метод расчёта объёма")
        
        # КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: Убедимся, что переменная method_var существует
        if not hasattr(self, 'method_var'):
            self.method_var = tk.StringVar(value="spline")
        
        methods = [
            ("Интеграл сплайна (1001 точка, рекоменд.)", "spline"),
            ("Метод Симпсона (501 точка)", "simpson"),
            ("Метод трапеций (2000 точек)", "trapezoidal"),
            ("Метод дисков (по точкам)", "disks"),
            ("Метод усечённых конусов", "frustums"),
        ]
        
        for text, value in methods:
            rb = ttk.Radiobutton(method_card, text=text, variable=self.method_var,
                               value=value, command=self.on_method_change)
            rb.pack(anchor='w', pady=2, padx=5)
        
        # Карточка: Уровень заполнения
        level_card = self.create_card(parent, "📊 Уровень заполнения")
        
        y_frame = ttk.Frame(level_card)
        y_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(y_frame, text="Уровень (см):").pack(side=tk.LEFT)
        self.y_level_var = tk.DoubleVar(value=0.0)
        y_entry = ttk.Entry(y_frame, textvariable=self.y_level_var, width=12)
        y_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(y_frame, text="Применить", 
                  command=self.apply_y_level).pack(side=tk.LEFT)
        
        self.y_slider = ttk.Scale(level_card, from_=0, to=100, 
                                 orient=tk.HORIZONTAL, length=300)
        self.y_slider.pack(fill=tk.X, pady=5)
        self.y_slider.bind('<ButtonRelease-1>', self.on_y_slider_release)
        
        percent_frame = ttk.Frame(level_card)
        percent_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(percent_frame, text="Заполнение (%):").pack(side=tk.LEFT)
        self.percent_var = tk.DoubleVar(value=0.0)
        percent_entry = ttk.Entry(percent_frame, textvariable=self.percent_var, width=12)
        percent_entry.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(percent_frame, text="Применить", 
                  command=self.apply_percent).pack(side=tk.LEFT)
        
        # Карточка: Результаты
        results_card = self.create_card(parent, "📈 Результаты расчёта")
        
        # Современное отображение результатов
        self.results_container = ttk.Frame(results_card)
        self.results_container.pack(fill=tk.BOTH, expand=True)
        
        # Создаем современный вид результатов
        self.create_modern_results_display(self.results_container)
        
        # Карточка: Быстрые действия
        actions_card = self.create_card(parent, "⚡ Быстрые действия")
        
        actions = [
            ("⚖️ Сравнить методы", self.compare_all_methods),
            ("💾 Сохранить профиль", self.save_current_profile),
            ("📈 Создать график", self.create_volume_chart),
            ("📋 Копировать результаты", self.copy_results_to_clipboard),
            ("🧪 Тест точности", self.test_volume_calculation)  # Новая кнопка
        ]
        
        for text, command in actions:
            btn = ttk.Button(actions_card, text=text, command=command)
            btn.pack(fill=tk.X, pady=3)
    
    def create_modern_results_display(self, parent):
        """Создание современного отображения результатов с КЛАССИФИКАЦИЕЙ ЦЕТЛИНА"""
        # Очищаем контейнер
        for widget in parent.winfo_children():
            widget.destroy()
        
        if not self.current_profile:
            label = ttk.Label(parent, text="Нет данных для отображения",
                            font=('Segoe UI', 11, 'italic'),
                            foreground=MODERN_PALETTE['secondary'])
            label.pack(expand=True)
            return
        
        # Основная информация
        main_frame = ttk.Frame(parent)
        main_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(main_frame, text=f"Профиль: {self.current_profile['name']}",
                 font=('Segoe UI', 10, 'bold'),
                 foreground=MODERN_PALETTE['primary']).pack(anchor='w')
        
        # Разделитель
        ttk.Separator(parent, orient='horizontal').pack(fill=tk.X, pady=5)
        
        # Метод расчета
        method_frame = ttk.Frame(parent)
        method_frame.pack(fill=tk.X, pady=5)
        
        method_names = {
            'spline': 'Интеграл сплайна (1001 точка)',
            'disks': 'Метод дисков',
            'frustums': 'Метод усечённых конусов',
            'trapezoidal': 'Метод трапеций (2000 точек)',
            'simpson': 'Метод Симпсона (501 точка)'
        }
        
        current_method = method_names.get(self.method_var.get(), self.method_var.get())
        ttk.Label(method_frame, text=f"Метод: {current_method}",
                 font=('Segoe UI', 9),
                 foreground=MODERN_PALETTE['dark']).pack(anchor='w')
        
        # Разделитель
        ttk.Separator(parent, orient='horizontal').pack(fill=tk.X, pady=5)
        
        # Объемы
        volumes_frame = ttk.Frame(parent)
        volumes_frame.pack(fill=tk.X, pady=5)
        
        # Заголовок объемов
        ttk.Label(volumes_frame, text="ОБЪЁМЫ",
                 font=('Segoe UI', 10, 'bold'),
                 foreground=MODERN_PALETTE['primary_dark']).pack(anchor='w', pady=(0, 5))
        
        if self.volume_calculator:
            try:
                # Получаем текущий уровень
                level = self.y_level_var.get()
                
                # Используем текущий метод расчета
                method = self.method_var.get()
                
                print(f"DEBUG: Расчет объема методом '{method}' до уровня {level} см")
                
                # Вычисляем объемы РАЗНЫМИ МЕТОДАМИ
                full_volume = self.volume_calculator.calculate_volume(method)
                level_volume = self.volume_calculator.calculate_volume(method, level)
                
                # Вычисляем процент заполнения
                if full_volume > 0:
                    percent = (level_volume / full_volume * 100)
                else:
                    percent = 0.0
                
                # Обновляем объем в профиле
                self.current_profile['volume'] = full_volume
                
                # Полный объем
                full_frame = ttk.Frame(volumes_frame)
                full_frame.pack(fill=tk.X, pady=2)
                
                ttk.Label(full_frame, text="🎯 Полный объём:",
                         font=('Segoe UI', 9, 'bold')).pack(side=tk.LEFT)
                ttk.Label(full_frame, text=f"{full_volume/1000:.3f} л ({full_volume:.1f} см³)",
                         font=('Segoe UI', 9)).pack(side=tk.RIGHT)
                
                # Объем до уровня
                level_frame = ttk.Frame(volumes_frame)
                level_frame.pack(fill=tk.X, pady=2)
                
                ttk.Label(level_frame, text=f"📏 До уровня {level:.1f} см:",
                         font=('Segoe UI', 9, 'bold')).pack(side=tk.LEFT)
                ttk.Label(level_frame, text=f"{level_volume/1000:.3f} л ({level_volume:.1f} см³)",
                         font=('Segoe UI', 9)).pack(side=tk.RIGHT)
                
                # Заполнение
                fill_frame = ttk.Frame(volumes_frame)
                fill_frame.pack(fill=tk.X, pady=2)
                
                ttk.Label(fill_frame, text="📈 Заполнение:",
                         font=('Segoe UI', 9, 'bold')).pack(side=tk.LEFT)
                ttk.Label(fill_frame, text=f"{percent:.1f}%",
                         font=('Segoe UI', 9)).pack(side=tk.RIGHT)
                
                # Геометрия
                geometry_frame = ttk.Frame(volumes_frame)
                geometry_frame.pack(fill=tk.X, pady=(10, 0))
                
                ttk.Label(geometry_frame, text="📐 Геометрия:",
                         font=('Segoe UI', 9, 'bold')).pack(anchor='w')
                
                height = np.max(self.current_profile['y'])
                diameter = np.max(self.current_profile['r']) * 2
                
                ttk.Label(geometry_frame, text=f"Высота: {height:.1f} см",
                         font=('Segoe UI', 9)).pack(anchor='w', padx=10)
                ttk.Label(geometry_frame, text=f"Макс. диаметр: {diameter:.1f} см",
                         font=('Segoe UI', 9)).pack(anchor='w', padx=10)
                
                # Классификация Цетлина
                if 'volume' in self.current_profile:
                    # Обновляем классификацию
                    tsetlin_classification = self.get_tsetlin_classification(self.current_profile['volume'])
                    self.current_profile['tsetlin_classification'] = tsetlin_classification
                    
                    tsetlin_info = tsetlin_classification
                    tsetlin_frame = ttk.Frame(volumes_frame)
                    tsetlin_frame.pack(fill=tk.X, pady=(10, 0))
                    
                    ttk.Label(tsetlin_frame, text="🎯 Классификация Цетлина:",
                             font=('Segoe UI', 9, 'bold')).pack(anchor='w')
                    
                    quality_color = 'green' if tsetlin_info['is_strict_quality'] else 'orange'
                    quality_text = "строгое" if tsetlin_info['is_strict_quality'] else "переходное"
                    
                    ttk.Label(tsetlin_frame, 
                             text=f"Группа {tsetlin_info['group']} ({tsetlin_info['group_name']})",
                             font=('Segoe UI', 9, 'bold'),
                             foreground=quality_color).pack(anchor='w', padx=10)
                    ttk.Label(tsetlin_frame, 
                             text=f"Качество: {quality_text}, Мобильность: {tsetlin_info['mobility_class']}",
                             font=('Segoe UI', 8)).pack(anchor='w', padx=20)
                
            except Exception as e:
                error_frame = ttk.Frame(volumes_frame)
                error_frame.pack(fill=tk.X, pady=5)
                
                ttk.Label(error_frame, text=f"⚠️ Ошибка расчета: {str(e)[:50]}...",
                         font=('Segoe UI', 9),
                         foreground=MODERN_PALETTE['danger']).pack()
        
        # Время расчета
        time_frame = ttk.Frame(parent)
        time_frame.pack(fill=tk.X, pady=(10, 0))
        
        ttk.Label(time_frame, text=f"✅ {datetime.now().strftime('%H:%M:%S')}",
                 font=('Segoe UI', 8),
                 foreground=MODERN_PALETTE['secondary']).pack(anchor='e')
    
    def create_card(self, parent, title):
        card = ttk.Frame(parent, style='Card.TFrame')
        card.pack(fill=tk.X, pady=8, padx=5)
        
        title_label = ttk.Label(card, text=title, style='Card.TLabel')
        title_label.pack(anchor='w', padx=10, pady=(10, 5))
        
        content = ttk.Frame(card)
        content.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0, 10))
        
        return content
    
    def create_status_bar(self, parent):
        self.status_var = tk.StringVar(value="Готов к работе")
        if HAVE_DND:
            self.status_var.set("Готов к работе. Перетащите DXF файлы из проводника!")
        
        status_bar = ttk.Label(parent, textvariable=self.status_var,
                              relief=tk.SUNKEN, anchor=tk.W,
                              background=MODERN_PALETTE['primary_light'],
                              foreground='white',
                              font=('Segoe UI', 10))
        status_bar.pack(fill=tk.X, pady=(10, 0))
    
    # ============================================================================
    # 3D УПРАВЛЕНИЕ МЕТОДЫ
    # ============================================================================
    
    def choose_3d_color(self):
        """Выбор цвета для 3D модели с мгновенным применением"""
        color = colorchooser.askcolor(title="Выберите цвет 3D модели", 
                                      initialcolor=self.surface_color_hex)
        if color[1]:
            self.surface_color_hex = color[1]
            self.color_button.config(bg=self.surface_color_hex)
            self.update_3d_plot()  # Мгновенное применение
    
    def reset_3d_view(self):
        """Сброс вида камеры к стандартному"""
        self.ax_3d.view_init(elev=30, azim=-60)
        self.ax_3d.set_proj_type(self.projection_type_3d_var.get())
        
        # Включаем/выключаем оси в зависимости от настройки
        if self.show_axes_3d_var.get():
            self.ax_3d.set_axis_on()
        else:
            self.ax_3d.set_axis_off()
        
        self.canvas_3d.draw()
    
    def set_isometric_view(self):
        """Установка изометрического вида"""
        self.ax_3d.view_init(elev=30, azim=45)
        self.ax_3d.set_proj_type(self.projection_type_3d_var.get())
        
        # Включаем/выключаем оси в зависимости от настройки
        if self.show_axes_3d_var.get():
            self.ax_3d.set_axis_on()
        else:
            self.ax_3d.set_axis_off()
        
        self.canvas_3d.draw()
    
    def set_top_view(self):
        """Установка вида сверху"""
        self.ax_3d.view_init(elev=90, azim=-90)
        self.ax_3d.set_proj_type(self.projection_type_3d_var.get())
        
        # Включаем/выключаем оси в зависимости от настройки
        if self.show_axes_3d_var.get():
            self.ax_3d.set_axis_on()
        else:
            self.ax_3d.set_axis_off()
        
        self.canvas_3d.draw()
    
    def save_3d_snapshot(self):
        """Сохранение снимка 3D модели"""
        if not self.current_profile:
            messagebox.showwarning("Ошибка", "Нет активной 3D модели для сохранения")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[
                ("PNG files", "*.png"),
                ("PDF files", "*.pdf"),
                ("SVG files", "*.svg"),
                ("All files", "*.*")
            ],
            title="Сохранить 3D модель",
            initialfile=f"{self.current_profile['name']}_3d.png"
        )
        
        if filename:
            try:
                self.fig_3d.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Успех", f"3D модель сохранена в {filename}")
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить: {str(e)}")
    
    def export_3d_model(self):
        """Экспорт 3D модели в формат STL"""
        if not HAVE_STL:
            messagebox.showerror("Ошибка экспорта", 
                "Для экспорта в STL необходима библиотека 'numpy-stl'. Установите её: pip install numpy-stl")
            return
        
        if self.X_surface is None or self.Y_surface is None or self.Z_surface is None:
            messagebox.showwarning("Ошибка", "Нет данных для экспорта. Постройте 3D модель сначала.")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".stl",
            filetypes=[("STL files", "*.stl"), ("All files", "*.*")],
            title="Сохранить 3D модель как STL",
            initialfile=f"{self.current_profile['name']}_3d.stl"
        )
        
        if filename:
            try:
                # Размеры сетки
                n_theta, n_y = self.X_surface.shape
                
                # Создаем список для хранения треугольников
                triangles = []
                
                for i in range(n_theta - 1):
                    for j in range(n_y - 1):
                        # Вершины четырехугольника
                        v1 = [self.X_surface[i, j], self.Y_surface[i, j], self.Z_surface[i, j]]
                        v2 = [self.X_surface[i+1, j], self.Y_surface[i+1, j], self.Z_surface[i+1, j]]
                        v3 = [self.X_surface[i+1, j+1], self.Y_surface[i+1, j+1], self.Z_surface[i+1, j+1]]
                        v4 = [self.X_surface[i, j+1], self.Y_surface[i, j+1], self.Z_surface[i, j+1]]
                        
                        # Разбиваем на два треугольника: v1,v2,v3 и v1,v3,v4
                        triangles.append([v1, v2, v3])
                        triangles.append([v1, v3, v4])
                
                # Преобразуем в массив numpy
                triangles_array = np.array(triangles)
                
                # Создаем STL mesh
                vessel_mesh = mesh.Mesh(np.zeros(triangles_array.shape[0], dtype=mesh.Mesh.dtype))
                for i, triangle in enumerate(triangles_array):
                    for j in range(3):
                        vessel_mesh.vectors[i][j] = triangle[j]
                
                # Сохраняем
                vessel_mesh.save(filename)
                messagebox.showinfo("Успех", f"3D модель сохранена в {filename}")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось экспортировать модель: {str(e)}")
                logging.error(f"Ошибка экспорта STL: {e}")
    
    def update_3d_plot(self, *args):
        """ИСПРАВЛЕННАЯ 3D модель с центрированием и управлением осями"""
        if not self.current_profile:
            return
        
        self.ax_3d.clear()
        
        y = self.current_profile['y']
        r = self.current_profile['r']
        
        # Получаем данные профиля
        y_prof = np.array(y)
        r_prof = np.array(r)
        
        # Нормализуем высоту
        y_prof = y_prof - np.min(y_prof)
        
        # Упрощение профиля для ускорения 3D-рендеринга
        original_points = np.column_stack([y_prof, r_prof])
        
        if self.settings['enable_3d_optimization'] and len(original_points) > 200:
            epsilon = self.settings['rdp_epsilon']
            simplified_points = self.simplify_profile_rdp(original_points, epsilon)
            
            y_prof = simplified_points[:, 0]
            r_prof = simplified_points[:, 1]
            
            reduction = (1 - len(y_prof)/len(original_points)) * 100
            print(f"Профиль упрощен: {len(original_points)} -> {len(y_prof)} точек ({reduction:.1f}% сокращение)")
        
        # Определяем количество сегментов по углу
        n_theta = self.segments_theta_var.get()
        
        # Создаём углы
        theta = np.linspace(0, 2 * np.pi, n_theta)
        
        # Создаём сетку
        theta_grid, h_grid = np.meshgrid(theta, y_prof, indexing='ij')
        r_expanded = r_prof[np.newaxis, :]  # Размер (1, n_y)
        r_expanded = np.tile(r_expanded, (n_theta, 1))  # Размер (n_theta, n_y)
        
        # Преобразуем в декартовы координаты
        X = r_expanded * np.cos(theta_grid)
        Z = r_expanded * np.sin(theta_grid)
        Y = np.tile(y_prof, (n_theta, 1))  # Высота одинакова для всех углов
        
        # Сохраняем данные для экспорта
        self.X_surface = X
        self.Y_surface = Y  # Высота
        self.Z_surface = Z
        
        # Настраиваем тип проекции из переменной
        self.ax_3d.set_proj_type(self.projection_type_3d_var.get())
        
        # Включаем/выключаем оси в зависимости от настройки
        if self.show_axes_3d_var.get():
            self.ax_3d.set_axis_on()
        else:
            self.ax_3d.set_axis_off()
        
        # Получаем текущие настройки из переменных
        current_alpha = self.alpha_3d_var.get()
        current_style = self.surface_style_3d_var.get()
        current_density = self.density_var.get()
        
        # Ключевое исправление: правильный порядок осей
        if current_style == 'solid':
            # Рисуем поверхность с выбранным цветом
            rstride_val = max(1, int(len(y_prof) / 50 * current_density))
            cstride_val = max(1, int(n_theta / 30 * current_density))
            
            # Исправленный вызов: X, Z, Y
            surface = self.ax_3d.plot_surface(
                X, Z, Y,  # X, Z, Y вместо X, Y, Z
                color=self.surface_color_hex,
                rstride=rstride_val,
                cstride=cstride_val,
                alpha=current_alpha,
                linewidth=0.3,
                antialiased=True,
                shade=True
            )
                
        elif current_style == 'wireframe':
            # Рисуем каркасную модель
            rstride_val = max(1, int(len(y_prof) / 20))
            cstride_val = max(1, int(n_theta / 20))
            
            # Исправленный вызов: X, Z, Y
            self.ax_3d.plot_wireframe(
                X, Z, Y,  # X, Z, Y вместо X, Y, Z
                color=self.surface_color_hex,
                rstride=rstride_val,
                cstride=cstride_val,
                alpha=current_alpha,
                linewidth=0.8
            )
        
        # Настройка осей в соответствии с новой системой координат
        self.ax_3d.set_xlabel('X (см)', color=MODERN_PALETTE['primary_dark'], fontsize=10)
        self.ax_3d.set_ylabel('Z (см)', color=MODERN_PALETTE['primary_dark'], fontsize=10)
        self.ax_3d.set_zlabel('Высота (см)', color=MODERN_PALETTE['primary_dark'], fontsize=10)
        self.ax_3d.set_title('3D модель сосуда', 
                           fontsize=12, color=MODERN_PALETTE['primary'])
        
        # Устанавливаем равный масштаб для всех осей (центрирование исправлено)
        self.set_axes_equal(self.ax_3d)
        
        # Устанавливаем угол обзора по умолчанию
        self.ax_3d.view_init(elev=30, azim=-60)
        
        # Добавляем сетку для лучшей ориентации
        self.ax_3d.grid(True, alpha=0.3)
        
        # Устанавливаем соотношение сторон для правильного отображения
        self.ax_3d.set_box_aspect([1, 1, 1])
        
        self.canvas_3d.draw()
        
        # Обновляем информацию о модели
        self.update_model_info()

    def update_model_info(self):
        """Обновление информации о модели"""
        if self.current_profile and self.model_info_label:
            height = np.max(self.current_profile['y'])
            diameter = np.max(self.current_profile['r']) * 2
            
            # Добавляем информацию о классификации Цетлина, если есть
            tsetlin_text = ""
            if 'tsetlin_classification' in self.current_profile:
                tsetlin_info = self.current_profile['tsetlin_classification']
                tsetlin_text = f"\n🎯 Группа Цетлина: {tsetlin_info['group']} ({tsetlin_info['group_name']})"
            
            info_text = f"📐 Модель: {self.current_profile['name']}\n"
            info_text += f"📏 Высота: {height:.1f} см\n"
            info_text += f"📏 Диаметр: {diameter:.1f} см{tsetlin_text}\n"
            info_text += f"🎨 Цвет: {self.surface_color_hex}\n"
            info_text += f"🎯 Стиль: {self.surface_style_3d_var.get()}"
            self.model_info_label.config(text=info_text)
    
    def show_performance_settings(self):
        """Показать окно настроек производительности"""
        settings_window = tk.Toplevel(self.root)
        settings_window.title("Настройки производительности")
        settings_window.geometry("400x300")
        settings_window.transient(self.root)
        settings_window.grab_set()
        
        main_frame = ttk.Frame(settings_window, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(main_frame, text="⚙️ Настройки производительности", 
                 font=('Segoe UI', 14, 'bold'),
                 foreground=MODERN_PALETTE['primary']).pack(pady=(0, 20))
        
        # Параметр упрощения RDP
        rdp_frame = ttk.Frame(main_frame)
        rdp_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(rdp_frame, text="Упрощение профиля (ε):", 
                 width=20).pack(side=tk.LEFT)
        rdp_var = tk.DoubleVar(value=self.settings['rdp_epsilon'])
        rdp_scale = ttk.Scale(rdp_frame, from_=0.001, to=0.1, 
                            variable=rdp_var, orient=tk.HORIZONTAL)
        rdp_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=10)
        ttk.Label(rdp_frame, textvariable=rdp_var).pack(side=tk.LEFT)
        
        # Количество сегментов 3D
        segments_frame = ttk.Frame(main_frame)
        segments_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(segments_frame, text="Сегментов 3D:", 
                 width=20).pack(side=tk.LEFT)
        segments_var = tk.IntVar(value=self.settings['3d_segments'])
        segments_spin = ttk.Spinbox(segments_frame, from_=10, to=100, 
                                  textvariable=segments_var, width=10)
        segments_spin.pack(side=tk.LEFT, padx=10)
        
        # Оптимизация
        opt_frame = ttk.Frame(main_frame)
        opt_frame.pack(fill=tk.X, pady=5)
        
        opt_var = tk.BooleanVar(value=self.settings['enable_3d_optimization'])
        ttk.Checkbutton(opt_frame, text="Включить оптимизацию 3D",
                       variable=opt_var).pack(anchor='w')
        
        # Кнопки
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=tk.X, pady=20)
        
        def apply_settings():
            self.settings['rdp_epsilon'] = rdp_var.get()
            self.settings['3d_segments'] = segments_var.get()
            self.settings['enable_3d_optimization'] = opt_var.get()
            
            # Обновить 3D модель если есть текущий профиль
            if self.current_profile:
                self.update_3d_plot()
            
            messagebox.showinfo("Настройки", "Настройки применены")
            settings_window.destroy()
        
        ttk.Button(btn_frame, text="Применить", 
                  command=apply_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Отмена", 
                  command=settings_window.destroy).pack(side=tk.LEFT, padx=5)
    
    # ============================================================================
    # ОСНОВНЫЕ ФУНКЦИИ
    # ============================================================================
    
    def add_dxf_files(self):
        files = filedialog.askopenfilenames(
            title="Выберите DXF файлы",
            filetypes=[("DXF files", "*.dxf"), ("All files", "*.*")]
        )
        
        if files:
            self.add_files_to_current_group(files)
    
    def create_new_group(self):
        dialog = tk.Toplevel(self.root)
        dialog.title("Новая группа")
        dialog.geometry("300x150")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="Введите название группы:").pack(pady=10)
        
        name_var = tk.StringVar()
        entry = ttk.Entry(dialog, textvariable=name_var, width=30)
        entry.pack(pady=5)
        entry.focus_set()
        
        def create():
            name = name_var.get().strip()
            if name and name not in self.groups:
                self.groups[name] = ProfileGroup(name)
                self.current_group = name
                self.update_tree()
                dialog.destroy()
        
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=10)
        
        ttk.Button(btn_frame, text="Создать", command=create).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Отмена", command=dialog.destroy).pack(side=tk.LEFT, padx=5)
    
    def rename_group(self):
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        if 'group' not in item.get('tags', []):
            return
        
        old_name = item['text']
        
        dialog = tk.Toplevel(self.root)
        dialog.title("Переименовать группу")
        dialog.geometry("300x150")
        
        ttk.Label(dialog, text="Новое название группы:").pack(pady=10)
        
        name_var = tk.StringVar(value=old_name)
        entry = ttk.Entry(dialog, textvariable=name_var, width=30)
        entry.pack(pady=5)
        entry.focus_set()
        entry.select_range(0, tk.END)
        
        def rename():
            new_name = name_var.get().strip()
            if new_name and new_name != old_name:
                if new_name in self.groups:
                    messagebox.showerror("Ошибка", f"Группа '{new_name}' уже существует")
                else:
                    group = self.groups[old_name]
                    group.name = new_name
                    self.groups[new_name] = group
                    del self.groups[old_name]
                    
                    if old_name in self.expanded_groups:
                        self.expanded_groups.remove(old_name)
                        self.expanded_groups.add(new_name)
                    
                    self.update_tree()
                    dialog.destroy()
        
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=10)
        
        ttk.Button(btn_frame, text="Переименовать", command=rename).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Отмена", command=dialog.destroy).pack(side=tk.LEFT, padx=5)
    
    def delete_group(self):
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        if 'group' not in item.get('tags', []):
            return
        
        group_name = item['text']
        
        if messagebox.askyesno("Подтверждение", f"Удалить группу '{group_name}' со всеми файлами?"):
            group = self.groups[group_name]
            
            # Удаляем файлы из общего списка
            for file_path in group.files:
                if file_path in self.profiles:
                    del self.profiles[file_path]
            
            # Удаляем группу
            del self.groups[group_name]
            
            if group_name in self.expanded_groups:
                self.expanded_groups.remove(group_name)
            
            if group_name == self.current_group:
                self.current_group = None
            
            self.update_tree()
            self.update_results_table()
            self.update_results_charts()
            
            if self.current_profile and self.current_profile['file_path'] not in self.profiles:
                self.current_profile = None
                self.volume_calculator = None
                self.update_profile_plot()
                self.update_3d_plot()
                self.create_modern_results_display(self.results_container)
    
    def sort_groups_by_name(self):
        sorted_groups = dict(sorted(self.groups.items()))
        self.groups = sorted_groups
        self.update_tree()
    
    def on_tree_select(self, event):
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        tags = item.get('tags', [])
        
        if 'file' in tags:
            file_path = item['values'][0]
            if file_path in self.profiles:
                self.display_profile(file_path)
    
    def on_tree_double_click(self, event):
        item = self.tree.identify_row(event.y)
        if item:
            item_data = self.tree.item(item)
            if 'group' in item_data.get('tags', []):
                group_name = item_data['text']
                if group_name in self.expanded_groups:
                    self.expanded_groups.remove(group_name)
                    self.tree.item(item, open=False)
                else:
                    self.expanded_groups.add(group_name)
                    self.tree.item(item, open=True)
    
    def show_tree_menu(self, event):
        item = self.tree.identify_row(event.y)
        if item:
            item_data = self.tree.item(item)
            tags = item_data.get('tags', [])
            
            if 'group' in tags:
                self.group_menu.post(event.x_root, event.y_root)
            elif 'file' in tags:
                self.tree_menu.post(event.x_root, event.y_root)
            else:
                self.tree_menu.post(event.x_root, event.y_root)
        else:
            self.tree_menu.post(event.x_root, event.y_root)
    
    def move_to_group(self):
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        if 'file' not in item.get('tags', []):
            return
        
        file_path = item['values'][0]
        
        dialog = tk.Toplevel(self.root)
        dialog.title("Переместить в группу")
        dialog.geometry("300x150")
        
        ttk.Label(dialog, text="Выберите группу:").pack(pady=10)
        
        group_var = tk.StringVar()
        groups = list(self.groups.keys())
        combo = ttk.Combobox(dialog, textvariable=group_var, values=groups, state='readonly')
        combo.pack(pady=5)
        
        def move():
            target_group = group_var.get()
            if target_group in self.groups:
                # Находим текущую группу файла
                source_group = None
                for group_name, group in self.groups.items():
                    if file_path in group.files:
                        source_group = group_name
                        break
                
                if source_group and source_group != target_group:
                    # Удаляем из старой группы
                    profile = self.groups[source_group].remove_profile(file_path)
                    # Добавляем в новую группу
                    self.groups[target_group].add_profile(profile, file_path)
                    
                    self.update_tree()
                    self.update_results_table()
                    self.update_results_charts()
                    dialog.destroy()
        
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(pady=10)
        
        ttk.Button(btn_frame, text="Переместить", command=move).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Отмена", command=dialog.destroy).pack(side=tk.LEFT, padx=5)
    
    def delete_selected(self):
        selection = self.tree.selection()
        if not selection:
            return
        
        item = self.tree.item(selection[0])
        tags = item.get('tags', [])
        
        if 'file' in tags:
            file_path = item['values'][0]
            
            if messagebox.askyesno("Подтверждение", f"Удалить профиль '{os.path.basename(file_path)}'?"):
                # Удаляем из группы
                for group in self.groups.values():
                    group.remove_profile(file_path)
                
                # Удаляем из общего списка
                if file_path in self.profiles:
                    del self.profiles[file_path]
                
                # Если удаляемый файл - текущий профиль, сбрасываем
                if self.current_profile and self.current_profile['file_path'] == file_path:
                    self.current_profile = None
                    self.volume_calculator = None
                    self.update_profile_plot()
                    self.update_3d_plot()
                    self.create_modern_results_display(self.results_container)
                
                self.update_tree()
                self.update_results_table()
                self.update_results_charts()
    
    def update_tree(self):
        expanded = []
        for group_name, group in self.groups.items():
            if group_name in self.expanded_groups:
                expanded.append(group_name)
        
        for item in self.tree.get_children():
            self.tree.delete(item)
        
        for group_name, group in self.groups.items():
            group_id = self.tree.insert('', 'end', text=group_name, tags=('group',))
            
            for file_path in group.files:
                filename = os.path.basename(file_path)
                profile = self.profiles.get(file_path)
                
                if profile:
                    # ВАЖНОЕ ИСПРАВЛЕНИЕ: Рассчитываем объем текущим методом для отображения в дереве
                    calculator = CorrectVolumeCalculator(profile['y'], profile['r'])
                    volume = calculator.calculate_volume(self.method_var.get())
                    
                    height = np.max(profile.get('y', [0]))
                    
                    # Добавляем информацию о классификации Цетлина, если есть
                    tsetlin_text = ""
                    if 'tsetlin_classification' in profile:
                        tsetlin_info = profile['tsetlin_classification']
                        tsetlin_text = f", Гр.{tsetlin_info['group']}"
                    
                    text = f"{filename} ({volume/1000:.2f} л{tsetlin_text}, H={height:.1f} см)"
                else:
                    text = f"{filename} (не обработан)"
                
                self.tree.insert(group_id, 'end', text=text, 
                               values=(file_path,), tags=('file',))
        
        for group_name in expanded:
            for item in self.tree.get_children():
                if self.tree.item(item, 'text') == group_name:
                    self.tree.item(item, open=True)
                    break
    
    def update_results_table(self):
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        
        for file_path, profile in self.profiles.items():
            if profile:
                group_name = self.find_profile_group(file_path)
                method = self.method_var.get()  # Используем текущий метод
                
                calculator = CorrectVolumeCalculator(profile['y'], profile['r'])
                volume = calculator.calculate_volume(method)  # Используем выбранный метод
                
                height = np.max(profile['y'])
                diameter = np.max(profile['r']) * 2
                
                method_names = {
                    'spline': 'Сплайн',
                    'disks': 'Диски',
                    'frustums': 'Конусы',
                    'trapezoidal': 'Трапеции',
                    'simpson': 'Симпсон'
                }
                method_display = method_names.get(method, method)
                
                # Получаем классификацию Цетлина
                tsetlin_group = ""
                if 'tsetlin_classification' in profile:
                    tsetlin_info = profile['tsetlin_classification']
                    tsetlin_group = f"{tsetlin_info['group']} ({tsetlin_info['group_name']})"
                
                self.results_tree.insert('', 'end', values=(
                    profile['name'],
                    group_name,
                    f'{volume/1000:.3f}',
                    f'{volume:.1f}',
                    f'{height:.1f}',
                    f'{diameter:.1f}',
                    method_display,
                    tsetlin_group
                ))
    
    def update_results_charts(self):
        """УЛУЧШЕННЫЙ МЕТОД: Обновление графиков с цветами по Цетлину"""
        if not self.profiles:
            return
        
        self.ax_chart1.clear()
        self.ax_chart2.clear()
        self.ax_chart3.clear()
        self.ax_chart4.clear()
        
        profile_names = []
        volumes = []
        heights = []
        diameters = []
        height_diameter_ratios = []  # Соотношение высота/диаметр
        tsetlin_groups = []
        tsetlin_colors = []
        
        method = self.method_var.get()
        
        # Собираем данные для всех профилей
        for profile in self.profiles.values():
            if profile:
                profile_names.append(profile['name'][:15])
                calculator = CorrectVolumeCalculator(profile['y'], profile['r'])
                volume = calculator.calculate_volume(method)
                volumes.append(volume / 1000)  # в литры
                height = np.max(profile['y'])
                diameter = np.max(profile['r']) * 2
                heights.append(height)
                diameters.append(diameter)
                height_diameter_ratios.append(height/diameter if diameter > 0 else 0)
                
                # Получаем группу Цетлина и цвет
                if 'tsetlin_classification' in profile:
                    group = profile['tsetlin_classification']['group']
                    tsetlin_groups.append(group)
                    group_num = self.roman_to_int(group)
                    
                    # Градиент синего: группы I-VII - светлые, VIII-XIV - средние, XV-XX - темные
                    if group_num <= 7:
                        color = plt.cm.Blues(0.3 + (group_num-1)/20)
                    elif group_num <= 14:
                        color = plt.cm.Blues(0.5 + (group_num-8)/20)
                    else:
                        color = plt.cm.Blues(0.7 + (group_num-15)/20)
                    tsetlin_colors.append(color)
                else:
                    tsetlin_groups.append("N/A")
                    tsetlin_colors.append('#95a5a6')  # Серый для без группы
        
        # 1. ГРАФИК: Распределение объёмов (с цветами по Цетлину)
        if volumes:
            bars = self.ax_chart1.bar(range(len(volumes)), volumes, 
                                     color=tsetlin_colors, edgecolor='white', linewidth=0.5)
            self.ax_chart1.set_title('Распределение объёмов по Цетлину', 
                                   fontsize=12, color=MODERN_PALETTE['primary'])
            self.ax_chart1.set_ylabel('Объём (литры)', fontsize=10)
            self.ax_chart1.set_xlabel('Профили', fontsize=10)
            self.ax_chart1.set_xticks(range(len(profile_names)))
            self.ax_chart1.set_xticklabels(profile_names, rotation=45, ha='right', fontsize=8)
            self.ax_chart1.grid(True, alpha=0.3, axis='y')
            
            # Добавляем значения
            for bar, volume in zip(bars, volumes):
                height = bar.get_height()
                if height > 0:  # Не добавляем подпись для нулевых значений
                    self.ax_chart1.text(bar.get_x() + bar.get_width()/2., height + max(volumes)*0.01,
                                       f'{volume:.2f}', ha='center', va='bottom', fontsize=7, rotation=90)
        
        # 2. ГРАФИК: Высота vs Диаметр (с цветами по Цетлину)
        if heights and diameters:
            scatter = self.ax_chart2.scatter(diameters, heights, 
                                            c=tsetlin_colors, s=60, alpha=0.7,
                                            edgecolors='white', linewidth=0.5)
            self.ax_chart2.set_title('Высота vs Диаметр по группам Цетлина', 
                                   fontsize=12, color=MODERN_PALETTE['primary'])
            self.ax_chart2.set_ylabel('Высота (см)', fontsize=10)
            self.ax_chart2.set_xlabel('Диаметр (см)', fontsize=10)
            self.ax_chart2.grid(True, alpha=0.3)
            
            # Линия тренда
            if len(diameters) > 1:
                z = np.polyfit(diameters, heights, 1)
                p = np.poly1d(z)
                x_trend = np.linspace(min(diameters), max(diameters), 100)
                self.ax_chart2.plot(x_trend, p(x_trend), "r--", alpha=0.5, linewidth=1, label='Линия тренда')
                self.ax_chart2.legend(fontsize=8)
        
        # 3. ГРАФИК: Соотношение высота/диаметр (с цветами по Цетлину)
        if height_diameter_ratios and tsetlin_groups:
            # Создаем список кортежей для сортировки
            data = list(zip(profile_names, height_diameter_ratios, tsetlin_groups, tsetlin_colors))
            # Сортируем по соотношению
            data.sort(key=lambda x: x[1])
            
            # Разделяем отсортированные данные
            sorted_names = [x[0] for x in data]
            sorted_ratios = [x[1] for x in data]
            sorted_groups = [x[2] for x in data]
            sorted_colors = [x[3] for x in data]
            
            bars = self.ax_chart3.bar(range(len(data)), sorted_ratios, 
                                     color=sorted_colors, edgecolor='white', linewidth=0.5)
            self.ax_chart3.set_title('Соотношение Высота/Диаметр по группам Цетлина', 
                                   fontsize=12, color=MODERN_PALETTE['primary'])
            self.ax_chart3.set_ylabel('Высота/Диаметр', fontsize=10)
            self.ax_chart3.set_xlabel('Профили', fontsize=10)
            self.ax_chart3.set_xticks(range(len(data)))
            self.ax_chart3.set_xticklabels(sorted_names, rotation=45, ha='right', fontsize=7)
            self.ax_chart3.grid(True, alpha=0.3, axis='y')
            
            # Средняя линия
            mean_ratio = np.mean(height_diameter_ratios)
            self.ax_chart3.axhline(y=mean_ratio, color='red', linestyle='--', alpha=0.7, 
                                  label=f'Среднее: {mean_ratio:.2f}')
            self.ax_chart3.legend(fontsize=8)
            
            # Подписи групп
            for i, (name, ratio, group, color) in enumerate(data):
                if group != 'N/A':
                    self.ax_chart3.text(i, ratio + 0.05, 
                                       f"Гр.{group}", 
                                       ha='center', va='bottom', fontsize=7, rotation=0)
        
        # 4. ГРАФИК: Классификация по Цетлину (только используемые группы)
        if tsetlin_groups and any(g != "N/A" for g in tsetlin_groups):
            group_counter = Counter([g for g in tsetlin_groups if g != "N/A"])
            
            if group_counter:
                # Сортируем группы по номеру
                sorted_groups = sorted(group_counter.keys(), key=lambda x: self.roman_to_int(x))
                group_counts = [group_counter[g] for g in sorted_groups]
                
                # Цвета для групп
                group_colors = []
                for group in sorted_groups:
                    group_num = self.roman_to_int(group)
                    if group_num <= 7:
                        color = plt.cm.Blues(0.3 + (group_num-1)/20)
                    elif group_num <= 14:
                        color = plt.cm.Blues(0.5 + (group_num-8)/20)
                    else:
                        color = plt.cm.Blues(0.7 + (group_num-15)/20)
                    group_colors.append(color)
                
                bars = self.ax_chart4.bar(range(len(sorted_groups)), group_counts,
                                         color=group_colors, edgecolor='white', linewidth=1)
                
                self.ax_chart4.set_title('Классификация сосудов по Цетлину', 
                                       fontsize=12, color=MODERN_PALETTE['primary'])
                self.ax_chart4.set_ylabel('Количество сосудов', fontsize=10)
                self.ax_chart4.set_xlabel('Группа качества', fontsize=10)
                self.ax_chart4.set_xticks(range(len(sorted_groups)))
                self.ax_chart4.set_xticklabels(sorted_groups, rotation=0)
                self.ax_chart4.grid(True, alpha=0.3, axis='y')
                
                # Значения на столбцах
                for bar, count in zip(bars, group_counts):
                    height = bar.get_height()
                    self.ax_chart4.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                                       str(count), ha='center', va='bottom', fontsize=9, fontweight='bold')
                
                # Статистика
                total_vessels = sum(group_counts)
                unique_groups = len(sorted_groups)
                self.ax_chart4.text(0.02, 0.98, 
                                   f'Всего: {total_vessels} сосудов\nГрупп: {unique_groups}',
                                   transform=self.ax_chart4.transAxes,
                                   fontsize=9, verticalalignment='top',
                                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        self.fig_charts.tight_layout()
        self.canvas_charts.draw()
    
    def process_files(self):
        unprocessed = []
        for file_path in self.profiles:
            if self.profiles[file_path] is None:
                unprocessed.append(file_path)
        
        if not unprocessed:
            messagebox.showinfo("Информация", "Все файлы уже обработаны")
            return
        
        thread = threading.Thread(target=self.process_files_thread, args=(unprocessed,))
        thread.daemon = True
        thread.start()
    
    def process_files_thread(self, files):
        for i, file_path in enumerate(files):
            self.status_var.set(f"Обработка: {os.path.basename(file_path)}...")
            
            try:
                profile = self.extract_profile_corrected(file_path)
                
                if profile:
                    self.profiles[file_path] = profile
                    
                    # Добавляем классификацию Цетлина
                    if 'volume' in profile:
                        tsetlin_classification = self.get_tsetlin_classification(profile['volume'])
                        profile['tsetlin_classification'] = tsetlin_classification
                    
                    self.root.after(0, self.update_tree)
                    self.root.after(0, self.update_results_table)
                    
                    if i == 0:
                        self.root.after(0, lambda: self.display_profile(file_path))
            
            except Exception as e:
                print(f"Ошибка обработки {file_path}: {e}")
        
        self.root.after(0, lambda: self.status_var.set("Обработка завершена"))
        self.root.after(0, self.update_results_charts)
    
    def extract_profile_corrected(self, file_path):
        try:
            doc = ezdxf.readfile(file_path)
            msp = doc.modelspace()
            
            points = []
            for entity in msp:
                if entity.dxftype() == 'LINE':
                    points.append(entity.dxf.start[:2])
                    points.append(entity.dxf.end[:2])
                elif entity.dxftype() in ['LWPOLYLINE', 'POLYLINE']:
                    try:
                        pts = entity.get_points()
                        for p in pts:
                            points.append(p[:2])
                    except:
                        pass
            
            if len(points) < 10:
                return None
            
            points = np.array(points)
            
            x_coords = points[:, 0]
            y_coords = points[:, 1]
            
            axis_x = np.min(x_coords)
            points[:, 0] -= axis_x
            
            y_min = np.min(points[:, 1])
            points[:, 1] -= y_min
            
            points *= 0.1
            
            radii = points[:, 0]
            heights = points[:, 1]
            
            if np.any(radii < -0.001):
                radii = np.abs(radii)
            
            sort_idx = np.argsort(heights)
            heights = heights[sort_idx]
            radii = radii[sort_idx]
            
            unique_heights, unique_idx = np.unique(heights, return_index=True)
            unique_radii = radii[unique_idx]
            
            if unique_heights[0] > 0.01:
                unique_heights = np.insert(unique_heights, 0, 0.0)
                unique_radii = np.insert(unique_radii, 0, unique_radii[0])
            
            n_points = 200
            if len(unique_heights) > 1:
                interp_func = interp1d(unique_heights, unique_radii, 
                                     kind='cubic', fill_value='extrapolate')
                
                max_height = np.max(unique_heights)
                interp_heights = np.linspace(0, max_height, n_points)
                interp_radii = interp_func(interp_heights)
                
                interp_radii = np.maximum(interp_radii, 0.0)
            else:
                interp_heights = np.array([0.0, 1.0])
                interp_radii = np.array([unique_radii[0], unique_radii[0]])
            
            # КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: Используем метод по умолчанию (диски) для начального расчета
            calculator = CorrectVolumeCalculator(interp_heights, interp_radii)
            
            # Используем метод дисков для начального расчета (позже пересчитается текущим методом)
            volume = calculator.method_disks()  # Исправлено: был method_spline_integral()
            
            # Создаем профиль
            profile = {
                'name': os.path.basename(file_path),
                'y': interp_heights,
                'r': interp_radii,
                'volume': volume,  # Сохраняем начальный объем
                'file_path': file_path,
                'is_half': True,
                'axis_x': axis_x
            }
            
            # Добавляем классификацию Цетлина
            tsetlin_classification = self.get_tsetlin_classification(volume)
            profile['tsetlin_classification'] = tsetlin_classification
            
            return profile
            
        except Exception as e:
            print(f"Ошибка извлечения профиля {file_path}: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def display_profile(self, file_path):
        profile = self.profiles.get(file_path)
        if not profile:
            return
        
        self.current_profile = profile
        self.volume_calculator = CorrectVolumeCalculator(profile['y'], profile['r'])
        
        self.update_profile_plot()
        self.update_3d_plot()
        self.update_volume_info()  # ВАЖНОЕ ИСПРАВЛЕНИЕ: добавляем вызов для пересчета объема
        self.create_modern_results_display(self.results_container)
        self.update_results_charts()
        
        # Обновление информации о модели
        if self.model_info_label:
            height = np.max(profile['y'])
            diameter = np.max(profile['r']) * 2
            volume = profile.get('volume', 0) / 1000
            
            # Добавляем информацию о классификации Цетлина
            tsetlin_text = ""
            if 'tsetlin_classification' in profile:
                tsetlin_info = profile['tsetlin_classification']
                tsetlin_text = f"\n🎯 Группа Цетлина: {tsetlin_info['group']} ({tsetlin_info['group_name']})"
            
            info_text = f"📐 Модель: {profile['name']}\n"
            info_text += f"📏 Высота: {height:.1f} см\n"
            info_text += f"📏 Диаметр: {diameter:.1f} см\n"
            info_text += f"🎯 Объём: {volume:.2f} л{tsetlin_text}\n"
            info_text += f"📊 Точек: {len(profile['y'])}"
            self.model_info_label.config(text=info_text)
        
        max_height = np.max(profile['y'])
        self.y_slider.config(to=max_height)
        self.y_level_var.set(0.0)
        self.y_slider.set(0.0)
    
    def update_profile_plot(self):
        if not self.current_profile:
            return
        
        self.ax_profile.clear()
        
        y = self.current_profile['y']
        r = self.current_profile['r']
        
        self.ax_profile.plot(r, y, color=MODERN_PALETTE['primary'], linewidth=2.5, label='Правая сторона')
        self.ax_profile.plot(-r, y, color=MODERN_PALETTE['primary'], linewidth=2.5, label='Левая сторона')
        self.ax_profile.axvline(x=0, color=MODERN_PALETTE['secondary'], linestyle='--', alpha=0.7, label='Ось симметрии')
        
        current_level = self.y_level_var.get()
        if current_level > 0:
            mask = y <= current_level
            
            self.ax_profile.fill_betweenx(y[mask], -r[mask], r[mask], 
                                        alpha=0.2, color=MODERN_PALETTE['primary_light'],
                                        label=f'Заполнение до {current_level:.1f} см')
            
            self.ax_profile.axhline(y=current_level, color=MODERN_PALETTE['accent'], 
                                 linestyle='-', linewidth=2, alpha=0.8)
            
            r_at_level = np.interp(current_level, y, r)
            self.ax_profile.plot(r_at_level, current_level, 'o', 
                               color=MODERN_PALETTE['accent'], markersize=8)
            self.ax_profile.plot(-r_at_level, current_level, 'o', 
                               color=MODERN_PALETTE['accent'], markersize=8)
        
        self.ax_profile.set_xlabel('Радиус (см)', fontsize=12, color=MODERN_PALETTE['primary_dark'])
        self.ax_profile.set_ylabel('Высота (см)', fontsize=12, color=MODERN_PALETTE['primary_dark'])
        
        # Добавляем информацию о классификации Цетлина в заголовок
        title_text = f"Профиль: {self.current_profile['name']}"
        if 'tsetlin_classification' in self.current_profile:
            tsetlin_info = self.current_profile['tsetlin_classification']
            title_text += f" | Группа Цетлина: {tsetlin_info['group']} ({tsetlin_info['group_name']})"
        
        self.ax_profile.set_title(title_text, 
                                fontsize=14, color=MODERN_PALETTE['primary'])
        self.ax_profile.legend(loc='upper right')
        self.ax_profile.grid(True, alpha=0.3, color='#D6DBDF')
        
        max_r = np.max(r)
        max_y = np.max(y)
        max_dim = max(max_r, max_y)
        
        self.ax_profile.set_xlim(-max_dim * 1.1, max_dim * 1.1)
        self.ax_profile.set_ylim(-max_y * 0.05, max_y * 1.05)
        self.ax_profile.set_aspect('equal', adjustable='box')
        
        self.canvas_profile.draw()
    
    def update_volume_info(self):
        if not self.volume_calculator or not self.current_profile:
            return
        
        method = self.method_var.get()
        level = self.y_level_var.get()
        
        try:
            # Используем правильный метод расчета
            full_volume = self.volume_calculator.calculate_volume(method)
            level_volume = self.volume_calculator.calculate_volume(method, level)
            
            # Защита от деления на ноль
            if full_volume > 0:
                percent = (level_volume / full_volume * 100)
            else:
                percent = 0.0
            
            self.percent_var.set(round(percent, 1))
            
            # ОБНОВЛЯЕМ объем в текущем профиле (убрали условие)
            self.current_profile['volume'] = full_volume
            
            # Обновляем классификацию Цетлина
            volume_cm3 = self.current_profile['volume']
            tsetlin_classification = self.get_tsetlin_classification(volume_cm3)
            self.current_profile['tsetlin_classification'] = tsetlin_classification
            
            # Обновляем отображение
            self.update_results_table()
            self.update_results_charts()
            self.create_modern_results_display(self.results_container)
            
        except Exception as e:
            print(f"Ошибка расчета объема методом {method}: {e}")
            # В случае ошибки используем метод дисков как резервный
            try:
                full_volume = self.volume_calculator.method_disks()
                level_volume = self.volume_calculator.method_disks(level)
                percent = (level_volume / full_volume * 100) if full_volume > 0 else 0
                
                self.percent_var.set(round(percent, 1))
                self.current_profile['volume'] = full_volume
                self.update_results_table()
                self.update_results_charts()
                self.create_modern_results_display(self.results_container)
            except Exception as e2:
                print(f"Резервный расчет тоже не удался: {e2}")
    
    def on_method_change(self):
        self.update_volume_info()
        self.update_profile_plot()
    
    def apply_y_level(self):
        try:
            level = float(self.y_level_var.get())
            max_height = np.max(self.current_profile['y']) if self.current_profile else 0
            
            if 0 <= level <= max_height:
                self.y_slider.set(level)
                self.update_volume_info()
                self.update_profile_plot()
            else:
                messagebox.showwarning("Ошибка", f"Уровень должен быть от 0 до {max_height:.1f}")
        except:
            messagebox.showwarning("Ошибка", "Введите корректное число")
    
    def on_y_slider_release(self, event):
        level = self.y_slider.get()
        self.y_level_var.set(round(level, 1))
        self.update_volume_info()
        self.update_profile_plot()
    
    def apply_percent(self):
        if not self.volume_calculator or not self.current_profile:
            return
        
        try:
            percent = float(self.percent_var.get())
            
            if 0 <= percent <= 100:
                method = self.method_var.get()
                
                full_volume = self.volume_calculator.calculate_volume(method)
                
                target_volume = full_volume * (percent / 100)
                
                max_height = np.max(self.current_profile['y'])
                low, high = 0, max_height
                mid = max_height / 2
                
                # Бинарный поиск уровня
                for _ in range(30):
                    mid = (low + high) / 2
                    mid_volume = self.volume_calculator.calculate_volume(method, mid)
                    
                    if abs(mid_volume - target_volume) < 0.1:
                        break
                    elif mid_volume < target_volume:
                        low = mid
                    else:
                        high = mid
                
                self.y_level_var.set(round(mid, 1))
                self.y_slider.set(mid)
                self.update_volume_info()
                self.update_profile_plot()
                
            else:
                messagebox.showwarning("Ошибка", "Процент должен быть от 0 до 100")
                
        except ValueError:
            messagebox.showwarning("Ошибка", "Введите корректное число")
    
    def on_profile_click(self, event):
        if event.inaxes != self.ax_profile:
            return
        
        level = event.ydata
        if level is not None:
            max_height = np.max(self.current_profile['y']) if self.current_profile else 0
            level = max(0, min(level, max_height))
            
            self.y_level_var.set(round(level, 1))
            self.y_slider.set(level)
            self.update_volume_info()
            self.update_profile_plot()
            
            self.dragging_level = True
    
    def on_profile_drag(self, event):
        if not self.dragging_level or event.inaxes != self.ax_profile:
            return
        
        level = event.ydata
        if level is not None:
            max_height = np.max(self.current_profile['y']) if self.current_profile else 0
            level = max(0, min(level, max_height))
            
            self.y_level_var.set(round(level, 1))
            self.y_slider.set(level)
            self.update_volume_info()
            self.update_profile_plot()
    
    def on_profile_release(self, event):
        self.dragging_level = False
    
    def compare_all_methods(self):
        if not self.volume_calculator:
            messagebox.showwarning("Ошибка", "Сначала загрузите профиль")
            return
        
        try:
            results = self.volume_calculator.calculate_all_methods()
            
            compare_window = tk.Toplevel(self.root)
            compare_window.title("Сравнение методов расчёта")
            compare_window.geometry("700x500")
            compare_window.configure(bg='white')
            
            tree_frame = ttk.Frame(compare_window)
            tree_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            
            tree = ttk.Treeview(tree_frame, columns=('Метод', 'Объём (л)', 'Объём (см³)', 'Отклонение'), 
                               show='headings', height=15)
            
            tree.heading('Метод', text='Метод')
            tree.heading('Объём (л)', text='Объём (л)')
            tree.heading('Объём (см³)', text='Объём (см³)')
            tree.heading('Отклонение', text='Отклонение от сплайна')
            
            tree.column('Метод', width=200)
            tree.column('Объём (л)', width=120)
            tree.column('Объём (см³)', width=120)
            tree.column('Отклонение', width=120)
            
            method_names = {
                'disks': 'Метод дисков',
                'frustums': 'Метод усечённых конусов',
                'trapezoidal': 'Метод трапеций',
                'simpson': 'Метод Симпсона',
                'spline': 'Интеграл сплайна (эталон)'
            }
            
            reference = results.get('spline', 0)
            
            for method, volume in results.items():
                if volume is not None and reference > 0:
                    deviation = ((volume - reference) / reference * 100)
                    
                    tree.insert('', 'end', values=(
                        method_names.get(method, method),
                        f'{volume/1000:.4f}',
                        f'{volume:.2f}',
                        f'{deviation:+.2f}%'
                    ))
            
            vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=tree.yview)
            tree.configure(yscrollcommand=vsb.set)
            
            tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            vsb.pack(side=tk.RIGHT, fill=tk.Y)
            
            # Добавляем информацию о точках
            info_frame = ttk.Frame(compare_window)
            info_frame.pack(fill=tk.X, padx=10, pady=(0, 10))
            
            info_text = f"""
            Информация о методах:
            • Сплайн: 1001 точка, интеграл Симпсона
            • Симпсон: 501 точка, интеграл Симпсона
            • Трапеции: 2000 точек, метод трапеций
            • Диски: интерполяция по точкам профиля
            • Конусы: усечённые конусы между точками
            """
            
            ttk.Label(info_frame, text=info_text, justify=tk.LEFT).pack()
            
            btn_frame = ttk.Frame(compare_window)
            btn_frame.pack(pady=10)
            
            ttk.Button(btn_frame, text="Закрыть", 
                      command=compare_window.destroy).pack()
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось сравнить методы: {str(e)}")
    
    def copy_results_to_clipboard(self):
        """Копирование результатов в буфер обмена с КЛАССИФИКАЦИЕЙ ЦЕТЛИНА"""
        if not self.current_profile:
            messagebox.showwarning("Ошибка", "Нет данных для копирования")
            return
        
        try:
            method = self.method_var.get()
            level = self.y_level_var.get()
            
            calculator = CorrectVolumeCalculator(self.current_profile['y'], self.current_profile['r'])
            full_volume = calculator.calculate_volume(method)
            level_volume = calculator.calculate_volume(method, level)
            percent = (level_volume / full_volume * 100) if full_volume > 0 else 0
            
            text = f"""Bobrinsky - Результаты анализа
Профиль: {self.current_profile['name']}
Метод расчета: {method}

Объемы:
- Полный объем: {full_volume/1000:.3f} л ({full_volume:.1f} см³)
- До уровня {level:.1f} см: {level_volume/1000:.3f} л ({level_volume:.1f} см³)
- Заполнение: {percent:.1f}%

Геометрия:
- Высота: {np.max(self.current_profile['y']):.1f} см
- Макс. диаметр: {np.max(self.current_profile['r']) * 2:.1f} см
"""
            
            # Добавляем классификацию Цетлина, если есть
            if 'tsetlin_classification' in self.current_profile:
                tsetlin_info = self.current_profile['tsetlin_classification']
                text += f"""
Классификация Цетлина:
- Группа качества: {tsetlin_info['group']} ({tsetlin_info['group_name']})
- Диапазон группы: {tsetlin_info['start_l']:.3f} – {tsetlin_info['end_l']:.3f} л
- Центр качества: {tsetlin_info['center_l']:.3f} л
- Класс мобильности: {tsetlin_info['mobility_class']}
- Строгое качество: {'Да' if tsetlin_info['is_strict_quality'] else 'Нет'}
"""
            
            text += f"""
Время расчета: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"""
            
            self.root.clipboard_clear()
            self.root.clipboard_append(text)
            messagebox.showinfo("Успех", "Результаты скопированы в буфер обмена")
            
        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось копировать: {str(e)}")
    
    def copy_tree_selection(self):
        selection = self.results_tree.selection()
        if not selection:
            return
        
        lines = []
        for item in selection:
            values = self.results_tree.item(item, 'values')
            lines.append('\t'.join(str(v) for v in values))
        
        text = '\n'.join(lines)
        self.root.clipboard_clear()
        self.root.clipboard_append(text)
    
    def show_tree_menu_results(self, event):
        self.tree_menu_results.post(event.x_root, event.y_root)
    
    def save_current_profile(self):
        if not self.current_profile:
            messagebox.showwarning("Ошибка", "Нет активного профиля")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile=f"{self.current_profile['name']}_profile.csv"
        )
        
        if filename:
            y = self.current_profile['y']
            r = self.current_profile['r']
            
            df = pd.DataFrame({
                'Высота_см': y,
                'Радиус_см': r,
                'Диаметр_см': r * 2,
                'Площадь_см2': np.pi * r**2
            })
            
            # Добавляем классификацию Цетлина в комментарии
            if 'tsetlin_classification' in self.current_profile:
                tsetlin_info = self.current_profile['tsetlin_classification']
                comments = [
                    f"# Классификация Цетлина для сосуда: {self.current_profile['name']}",
                    f"# Группа качества: {tsetlin_info['group']} ({tsetlin_info['group_name']})",
                    f"# Объем: {tsetlin_info['volume_l']:.3f} л",
                    f"# Диапазон группы: {tsetlin_info['start_l']:.3f} – {tsetlin_info['end_l']:.3f} л",
                    f"# Класс мобильности: {tsetlin_info['mobility_class']}",
                    f"# Строгое качество: {'Да' if tsetlin_info['is_strict_quality'] else 'Нет'}"
                ]
                
                with open(filename, 'w', encoding='utf-8') as f:
                    for comment in comments:
                        f.write(comment + '\n')
                    df.to_csv(f, index=False)
            else:
                df.to_csv(filename, index=False, encoding='utf-8')
            
            messagebox.showinfo("Успех", f"Профиль сохранён в {filename}")
    
    def create_volume_chart(self):
        if not self.volume_calculator:
            return
        
        chart_window = tk.Toplevel(self.root)
        chart_window.title("График сравнения методов")
        chart_window.geometry("800x600")
        
        fig = Figure(figsize=(10, 6), dpi=100)
        ax = fig.add_subplot(111)
        
        results = self.volume_calculator.calculate_all_methods()
        valid_results = {k: v for k, v in results.items() if v is not None}
        
        if not valid_results:
            return
        
        methods = list(valid_results.keys())
        volumes = list(valid_results.values())
        
        display_names = {
            'disks': 'Диски',
            'frustums': 'Конусы',
            'trapezoidal': 'Трапеции',
            'simpson': 'Симпсон',
            'spline': 'Сплайн'
        }
        
        display_methods = [display_names.get(m, m) for m in methods]
        
        bars = ax.bar(display_methods, [v/1000 for v in volumes], 
                     color=GRADIENT[:len(methods)])
        
        ax.set_title('Сравнение методов расчёта объёма', 
                   fontsize=14, color=MODERN_PALETTE['primary'])
        ax.set_ylabel('Объём (литры)', color=MODERN_PALETTE['primary_dark'])
        ax.set_xlabel('Метод расчёта', color=MODERN_PALETTE['primary_dark'])
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, volume in zip(bars, volumes):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                   f'{volume/1000:.3f} л', ha='center', va='bottom')
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, chart_window)
        canvas.draw()
        
        toolbar = CustomNavigationToolbar(canvas, chart_window)
        
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def export_excel(self):
        if not self.profiles:
            messagebox.showwarning("Ошибка", "Нет данных для экспорта")
            return
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                    profile_data = []
                    for file_path, profile in self.profiles.items():
                        if profile:
                            # КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: Рассчитываем объем текущим методом
                            calculator = CorrectVolumeCalculator(profile['y'], profile['r'])
                            volume = calculator.calculate_volume(self.method_var.get())
                            
                            # Получаем классификацию Цетлина
                            tsetlin_group = ""
                            tsetlin_name = ""
                            if 'tsetlin_classification' in profile:
                                tsetlin_info = profile['tsetlin_classification']
                                tsetlin_group = tsetlin_info['group']
                                tsetlin_name = tsetlin_info['group_name']
                            
                            profile_data.append({
                                'Имя файла': profile['name'],
                                'Группа': self.find_profile_group(file_path),
                                'Объём (л)': volume / 1000,
                                'Объём (см³)': volume,
                                'Высота (см)': np.max(profile.get('y', [0])),
                                'Диаметр (см)': np.max(profile.get('r', [0])) * 2,
                                'Точек': len(profile.get('y', [])),
                                'Группа Цетлина': tsetlin_group,
                                'Качество Цетлина': tsetlin_name,
                                'Метод расчёта': self.method_var.get()
                            })
                    
                    if profile_data:
                        df_profiles = pd.DataFrame(profile_data)
                        df_profiles.to_excel(writer, sheet_name='Профили', index=False)
                    
                    if self.current_profile:
                        y = self.current_profile['y']
                        r = self.current_profile['r']
                        
                        detail_data = {
                            'Высота_см': y,
                            'Радиус_см': r,
                            'Диаметр_см': r * 2,
                            'Площадь_см2': np.pi * r**2
                        }
                        
                        df_detail = pd.DataFrame(detail_data)
                        df_detail.to_excel(writer, sheet_name='Детали', index=False)
                    
                    # Добавляем лист с классификацией Цетлина
                    tsetlin_data = []
                    for class_data in self.tsetlin_classification:
                        tsetlin_data.append({
                            'Группа': class_data['group'],
                            'Начало (л)': class_data['start_l'],
                            'Центр (л)': class_data['center_l'],
                            'Конец (л)': class_data['end_l'],
                            'Качество': class_data['quality_name'],
                            'Класс мобильности': class_data['mobility_class'],
                            'Описание': class_data['description']
                        })
                    
                    df_tsetlin = pd.DataFrame(tsetlin_data)
                    df_tsetlin.to_excel(writer, sheet_name='Классификация Цетлина', index=False)
                
                messagebox.showinfo("Успех", f"Данные экспортированы в {filename}")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось экспортировать: {str(e)}")
    
    def find_profile_group(self, file_path):
        for group_name, group in self.groups.items():
            if file_path in group.files:
                return group_name
        return "Без группы"
    
    def update_plots(self):
        if self.current_profile:
            self.update_profile_plot()
            self.update_3d_plot()
            self.update_volume_info()
            self.update_results_charts()
    
    def show_tsetlin_info(self):
        """Показать подробную информацию о классификации Цетлина"""
        info_text = """🎯 НАУЧНАЯ КЛАССИФИКАЦИЯ СОСУДОВ ПО Ю.Б. ЦЕТЛИНУ

ОСНОВНЫЕ ПРИНЦИПЫ:
• Классификация основана на объеме сосуда (в литрах)
• Используется логарифмическая шкала с интервалами ±0.95
• Каждая качественная группа имеет центр, начало и конец
• Группы сгруппированы в 6 классов мобильности

КЛАССЫ МОБИЛЬНОСТИ:
1. Супермалые (< 0,097 л) - ароматические вещества
2. Мобильные (0,097–50,0 л) - легко перемещаются одним человеком
3. Ограниченно-мобильные (50,0–200,0 л) - требуют усилий 2+ человек
4. Мало-мобильные (200,0–800,0 л) - перемещались крайне редко
5. Условно-мобильные (800,0–3200,0 л) - только пустыми
6. Стационарные (3200,0–25000,0 л) - не предполагают перемещения

НАУЧНАЯ ЗНАЧИМОСТЬ:
• Позволяет определить функциональное назначение сосудов
• Отражает технологические возможности древних мастеров
• Дает представление о логистике и транспортировке
• Помогает в культурной и хронологической атрибуции

МЕТОДОЛОГИЯ:
1. Вычисляется точный объем сосуда (рекомендуется метод сплайна)
2. Объем сравнивается с эталонной шкалой Цетлина
3. Определяется качественная группа и класс мобильности
4. Анализируется положение объема относительно центра группы

В программе используется полная шкала из 20 качественных групп,
разработанная Ю.Б. Цетлиным для археологических исследований."""
        
        info_window = tk.Toplevel(self.root)
        info_window.title("Классификация Цетлина - Научная справка")
        info_window.geometry("600x500")
        info_window.configure(bg='white')
        
        text = tk.Text(info_window, wrap=tk.WORD, 
                      font=('Segoe UI', 11), 
                      bg='white',
                      fg=MODERN_PALETTE['dark'],
                      padx=20, pady=20)
        text.insert(1.0, info_text)
        text.config(state='disabled')
        
        scrollbar = ttk.Scrollbar(info_window, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        btn_frame = ttk.Frame(info_window)
        btn_frame.pack(fill=tk.X, pady=(0, 10), padx=20)
        
        ttk.Button(btn_frame, text="Закрыть", 
                  command=info_window.destroy).pack()
    
    def show_help(self):
        help_text = """🎯 BOBRINSKY - АНАЛИЗАТОР СОСУДОВ v6.0 (С КЛАССИФИКАЦИЕЙ ЦЕТЛИНА)

ОСНОВНЫЕ ВОЗМОЖНОСТИ:

1. 📊 НОВАЯ СТРУКТУРА ВКЛАДОК:
   • "Объём" - объединенная вкладка с подвкладками:
     - "Профиль и Объемы" - график профиля и расчеты
     - "Таблицы" - результаты анализа и шкала Цетлина
     - "Графики" - визуализация данных (с графиком Цетлина)
     - "Классификация Цетлина" - подробная информация

2. 🎯 КЛАССИФИКАЦИЯ ЦЕТЛИНА (НАУЧНАЯ):
   • Автоматическое определение качественной группы объема
   • 20 групп качества с центрами и интервалами
   • 6 классов мобильности с описаниями
   • Визуализация на графике и в таблицах
   • Экспорт классификации в Excel

3. 🏺 УЛУЧШЕННАЯ 3D-ВИЗУАЛИЗАЦИЯ:
   • ИСПРАВЛЕННОЕ ЦЕНТРИРОВАНИЕ - сосуд теперь в центре координат
   • Включение/выключение осей координат
   • Корректный масштаб по всем осям
   • Все настройки применяются мгновенно

4. 📏 ТОЧНЫЙ РАСЧЁТ ОБЪЁМА:
   • 5 методов расчёта с разным количеством точек
   • Интеграл сплайна (1001 точка) - рекомендованный метод
   • Управление уровнем заполнения
   • Сравнение методов

5. 📂 УПРАВЛЕНИЕ ДАННЫМИ:
   • Drag-and-drop DXF файлов (при установке tkinterdnd2)
   • Группировка профилей
   • Экспорт в Excel с классификацией Цетлина
   • Сохранение профилей в CSV

ИНСТРУКЦИЯ:
1. Добавьте DXF файлы (кнопка или drag-and-drop)
2. Обработайте файлы для извлечения профилей
3. Выберите профиль в дереве слева
4. Настройте уровень заполнения и метод расчёта
5. Перейдите на вкладку "Классификация Цетлина" для анализа
6. Используйте 3D-визуализацию для изучения геометрии
7. Экспортируйте результаты при необходимости

НАУЧНАЯ ТЕРМИНОЛОГИЯ:
• Качественная группа объема - вместо "категория размера"
• Класс мобильности - функциональная классификация
• Строгое качество - объем близок к центру группы
• Переходная зона - объем между группами

УСТАНОВКА БИБЛИОТЕК:
• Для drag-and-drop: pip install tkinterdnd2
• Для экспорта STL: pip install numpy-stl

Версия 6.0 включает полную реализацию научной классификации
сосудов по методике Ю.Б. Цетлина для археологических исследований."""
        
        help_window = tk.Toplevel(self.root)
        help_window.title("Справка - Bobrinsky v6.0")
        help_window.geometry("600x500")
        help_window.configure(bg='white')
        
        text = tk.Text(help_window, wrap=tk.WORD, 
                      font=('Segoe UI', 11), 
                      bg='white',
                      fg=MODERN_PALETTE['dark'],
                      padx=20, pady=20)
        text.insert(1.0, help_text)
        text.config(state='disabled')
        
        scrollbar = ttk.Scrollbar(help_window, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        btn_frame = ttk.Frame(help_window)
        btn_frame.pack(fill=tk.X, pady=(0, 10), padx=20)
        
        ttk.Button(btn_frame, text="Закрыть", 
                  command=help_window.destroy).pack()
    
    def test_volume_calculation(self):
        """Тестирование расчета объемов разными методами"""
        if not self.current_profile:
            messagebox.showwarning("Ошибка", "Сначала загрузите профиль")
            return
        
        print("\n" + "="*60)
        print("ТЕСТИРОВАНИЕ РАСЧЕТА ОБЪЕМОВ")
        print("="*60)
        
        y = self.current_profile['y']
        r = self.current_profile['r']
        
        print(f"Количество точек профиля: {len(y)}")
        print(f"Высота: от {np.min(y):.2f} до {np.max(y):.2f} см")
        print(f"Радиус: от {np.min(r):.2f} до {np.max(r):.2f} см")
        print(f"Максимальный диаметр: {np.max(r) * 2:.2f} см")
        
        calculator = CorrectVolumeCalculator(y, r)
        
        print("\n1. РАСЧЕТ ПОЛНОГО ОБЪЕМА ВСЕМИ МЕТОДАМИ:")
        print("-" * 50)
        results = calculator.calculate_all_methods()
        
        print("\n2. РАСЧЕТ ОБЪЕМА ДО РАЗНЫХ УРОВНЕЙ:")
        print("-" * 50)
        
        test_levels = [0.25, 0.5, 0.75]
        for level in test_levels:
            level_cm = level * np.max(y)
            print(f"\nУровень: {level_cm:.2f} см ({level*100:.0f}% высоты):")
            
            for method_name, method_desc in [
                ('disks', 'Метод дисков'),
                ('frustums', 'Метод усечённых конусов'),
                ('spline', 'Интеграл сплайна')
            ]:
                try:
                    vol = calculator.calculate_volume(method_name, level_cm)
                    print(f"  {method_desc}: {vol/1000:.6f} л")
                except Exception as e:
                    print(f"  {method_desc}: Ошибка - {e}")
        
        print("\n3. ТОЧНОСТЬ МЕТОДОВ:")
        print("-" * 50)
        
        if 'spline' in results and results['spline'] is not None:
            reference = results['spline']
            print(f"Эталонный объем (метод сплайна): {reference/1000:.6f} л")
            
            for method_name in ['disks', 'frustums', 'trapezoidal', 'simpson']:
                if method_name in results and results[method_name] is not None:
                    vol = results[method_name]
                    diff = vol - reference
                    diff_percent = (diff / reference) * 100
                    print(f"{method_name}: {vol/1000:.6f} л, разница: {diff_percent:+.3f}%")
        
        print("\n4. ВАЛИДАЦИЯ:")
        print("-" * 50)
        
        # Проверка на монотонность
        volumes_at_levels = []
        test_heights = np.linspace(0, np.max(y), 6)
        
        for h in test_heights:
            vol = calculator.calculate_volume('spline', h)
            volumes_at_levels.append(vol)
        
        is_monotonic = all(volumes_at_levels[i] <= volumes_at_levels[i+1] 
                          for i in range(len(volumes_at_levels)-1))
        
        print(f"Монотонность объемов: {'✓' if is_monotonic else '✗'}")
        print(f"Объем на 100% высоты: {volumes_at_levels[-1]/1000:.6f} л")
        
        # Сохранение результатов в файл
        try:
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"volume_test_{timestamp}.txt"
            
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(f"ТЕСТ РАСЧЕТА ОБЪЕМОВ - {timestamp}\n")
                f.write(f"Профиль: {self.current_profile['name']}\n")
                f.write(f"Высота: {np.max(y):.2f} см\n")
                f.write(f"Макс. диаметр: {np.max(r)*2:.2f} см\n\n")
                
                f.write("ПОЛНЫЕ ОБЪЕМЫ:\n")
                for method, volume in results.items():
                    if volume is not None:
                        f.write(f"{method}: {volume/1000:.6f} л\n")
                
                f.write("\nТОЧНОСТЬ МЕТОДОВ:\n")
                if 'spline' in results:
                    for method in ['disks', 'frustums', 'trapezoidal', 'simpson']:
                        if method in results and results[method] is not None:
                            diff = (results[method] - results['spline']) / results['spline'] * 100
                            f.write(f"{method}: {diff:+.3f}%\n")
            
            print(f"\nРезультаты сохранены в файл: {filename}")
            
        except Exception as e:
            print(f"Ошибка сохранения результатов: {e}")
    
    def lighten_color(self, color, factor=0.2):
        try:
            color = color.lstrip('#')
            r, g, b = int(color[0:2], 16), int(color[2:4], 16), int(color[4:6], 16)
            
            r = min(255, int(r + (255 - r) * factor))
            g = min(255, int(g + (255 - g) * factor))
            b = min(255, int(b + (255 - b) * factor))
            
            return f'#{r:02x}{g:02x}{b:02x}'
        except:
            return color
    
    def start_queue_processor(self):
        def process():
            try:
                while True:
                    task = self.processing_queue.get_nowait()
                    
                    if task[0] == 'update_status':
                        self.status_var.set(task[1])
                    
                    self.processing_queue.task_done()
            except queue.Empty:
                pass
            finally:
                self.root.after(100, process)
        
        self.root.after(100, process)

# ============================================================================
# ЗАПУСК ПРОГРАММЫ
# ============================================================================

def main():
    # Используем TkinterDnD если доступен, иначе стандартный tkinter
    if HAVE_DND:
        root = TkinterDnD.Tk()
    else:
        root = tk.Tk()
    
    app = BobrinskyAnalyzer(root)
    
    root.update_idletasks()
    width = root.winfo_width()
    height = root.winfo_height()
    x = (root.winfo_screenwidth() // 2) - (width // 2)
    y = (root.winfo_screenheight() // 2) - (height // 2)
    root.geometry(f'{width}x{height}+{x}+{y}')
    
    root.mainloop()

if __name__ == "__main__":
    main()