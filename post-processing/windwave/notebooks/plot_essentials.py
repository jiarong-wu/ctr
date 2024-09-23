
''' For shared colorbar '''

fig, axes = plt.subplots(1, 3, sharey=True, figsize=[5,2])

# Plotting ...

cbar_ax = fig.add_axes([1, 0.2, 0.01, 0.6])  # Adjust as needed
cbar = plt.colorbar(a, cax=cbar_ax)
cbar.set_label(r'$u_x$')

axes[0].set_ylabel('$y$', fontsize=6)
axes[0].set_yticks([0,1,2])

plt.tight_layout()