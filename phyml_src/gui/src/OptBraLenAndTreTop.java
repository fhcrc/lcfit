
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * Implements all components necessary for
 * "Tree topology search and "Optimise branch length".
 * 
 * @author Christoph Knapp
 */

public class OptBraLenAndTreTop extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox TreTopBox;
	private JComboBox OptBraLenBox;
	private JLabel lab1;
	private JLabel lab2;

	/**
	 * Implements all components necessary for
	 * "Tree topology search and "Optimise branch length".
	 * 
	 * @param isYes
	 *            Whether the DropDown menu for selecting Optimise Branch length
	 *            is visible (false) or not (true).
	 */
	public OptBraLenAndTreTop(boolean isYes) {
		lab1 = new JLabel("Tree Topology Search");
		lab2 = new JLabel("Opt. Branch Length");
		TreTopBox = new JComboBox(new String[] { "NNI", "SPR", "NNI and SPR" });
		OptBraLenBox = new JComboBox(new String[] { "yes", "no" });
		TreTopBox.addActionListener(this);
		if (isYes) {
			OptBraLenBox.setVisible(false);
			lab2.setVisible(false);
		} else {
			TreTopBox.setVisible(false);
			lab1.setVisible(false);
		}
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.44, 0.9);
		JPanel p1 = new JPanel();
		add(p1);
		layout.setDimensions(0.1, 0.9);
		add(new JPanel());
		layout.setDimensions(0.44, 0.9);
		JPanel p2 = new JPanel();
		add(p2);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.62, 1);
		p1.add(lab1);
		lO1.setDimensions(0.1, 1);
		p1.add(new JLabel());
		lO1.setDimensions(0.28, 1);
		p1.add(TreTopBox);
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.61, 1);
		p2.add(lab2);
		lO2.setDimensions(0.1, 1);
		p2.add(new JLabel());
		lO2.setDimensions(0.23, 1);
		p2.add(OptBraLenBox);
	}

	public void actionPerformed(ActionEvent e) {
		if (TreTopBox.getSelectedItem().toString().equals("NNI")) {
			PhymlPanel.tS.setIsNNI(true);
		} else {
			PhymlPanel.tS.setIsNNI(false);
		}
	}

	/**
	 * Sets whether the Optimise Branch length components or the Tree topology
	 * components are visible.
	 * 
	 * @param b
	 *            If true the optimise branch length components are not visible.
	 *            and the Tree topology components are visible.<br>
	 *            If false its the opposite from above.
	 */
	public void setIsYes(boolean b) {
		if (b) {
			OptBraLenBox.setVisible(false);
			lab2.setVisible(false);
			TreTopBox.setVisible(true);
			if (TreTopBox.getSelectedIndex() > 0) {
				PhymlPanel.tS.setIsNNI(false);
			} else {
				PhymlPanel.tS.setIsNNI(true);
			}
			lab1.setVisible(true);
		} else {
			OptBraLenBox.setVisible(true);
			lab2.setVisible(true);
			TreTopBox.setVisible(false);
			PhymlPanel.tS.setIsNNI(true);
			lab1.setVisible(false);
		}
	}

	/**
	 * Passes the user choice from the "Tree Topology Search" options.
	 * 
	 * @return "BEST" if "NNI and SPR" is selected and "SPR" if "SPR" is
	 *         selected, ""(empty String) otherwise.
	 */
	public String getSearch() {
		if (PhymlPanel.tS.isOptimiseTreeTopology()) {
			if (TreTopBox.getSelectedItem().toString().equals("NNI and SPR")) {
				return "BEST";
			} else if (TreTopBox.getSelectedItem().toString().equals("SPR")) {
				return TreTopBox.getSelectedItem().toString();
			}
		}
		return "";
	}

	/**
	 * Retrieves whether the "Optimise Branch Length" parameter is set to yes or
	 * no.
	 * 
	 * @return Returns true if "Optimise Branch Length" is set to yes otherwise
	 *         false.
	 */
	public boolean isOptimiseBranchLength() {
		if (OptBraLenBox.getSelectedItem().toString().equals("yes")) {
			return true;
		}
		return false;
	}

}
