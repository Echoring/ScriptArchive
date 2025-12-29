import matplotlib.pyplot as plt

labels_outer = ['0-1000', '1000-2500', 
                '2500-4000', '4000-5000']
sizes_outer = [31, 20, 17, 8]  # Total
labels_inner = ['ASEG', 'not ASEG'] * 4
sizes_inner = [24, 7, 14, 6, 11, 6, 3, 5]  # ASEG/not ASEG

fig, ax = plt.subplots(figsize=(10, 6))
ax.pie(sizes_outer, labels=labels_outer, radius=1.2, wedgeprops=dict(width=0.3, edgecolor='w'))
ax.pie(sizes_inner, labels=labels_inner, radius=0.9, wedgeprops=dict(width=0.3, edgecolor='w'))
ax.set(aspect="equal", title='Nested Pie Chart')
plt.show()

fig, ax = plt.subplots(figsize=(10, 6))
ax.pie(sizes_outer, labels=labels_outer, radius=1.2, wedgeprops=dict(width=0.3, edgecolor='w'))
ax.pie(sizes_inner, labels=labels_inner, radius=0.9, wedgeprops=dict(width=0.3, edgecolor='w'))
ax.set(aspect="equal", title='Nested Pie Chart')

fig.savefig('nested_pie2.svg', format='svg', bbox_inches='tight')