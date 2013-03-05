
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Class for implementing all components to set the values for whether
 * "Non parametric Bootstrapping" is used, how many "Replicates", whether
 * additional statistics are printed to file and which
 * "Approximate Likelihood Ratio Test" is used.
 * 
 * @author Christoph Knapp
 * 
 */
public class BooStrAndAprLikRatTes extends JPanel implements ActionListener,
		FocusListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox booStrBox;
	private JComboBox likRatBox;
	private JComboBox priBooBox;
	private CustomTextField numRepField;

	/**
	 * Constructor method for implementing all components to set the values for
	 * whether "Non parametric Bootstrapping" is used, how many "Replicates",
	 * whether additional statistics are printed to file and which
	 * "Approximate Likelihood Ratio Test" is used.
	 */
	public BooStrAndAprLikRatTes() {
		booStrBox = new JComboBox(new String[] { "no", "yes" });
		likRatBox = new JComboBox(new String[] { "no", "SH-like supports",
				"aLRT statistics", "Chi2-based supports" });
		priBooBox = new JComboBox(new String[] { "no", "yes" });
		likRatBox.setSelectedIndex(1);
		numRepField = new CustomTextField("");
		numRepField.setEnabled(false);
		numRepField.addFocusListener(this);
		priBooBox.setEnabled(false);
		booStrBox.addActionListener(this);
		likRatBox.addActionListener(this);
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.98, 0.36);
		JPanel p1 = new JPanel();
		add(p1);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.98, 0.1);
		add(new JPanel());
		layout.setDimensions(0.98, 0.36);
		JPanel p2 = new JPanel();
		add(p2);
		layout.setDimensions(98, 0.1);
		add(new JPanel());
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.32, 1);
		p1.add(new JLabel("Non Parametric Bootstrapping"));
		lO1.setDimensions(0.13, 1);
		p1.add(booStrBox);
		lO1.setDimensions(0.05, 1);
		p1.add(new JPanel());
		lO1.setDimensions(0.12, 1);
		p1.add(new JLabel("Replicates"));
		lO1.setDimensions(0.08, 1);
		p1.add(numRepField);
		lO1.setDimensions(0.05, 1);
		p1.add(new JPanel());
		lO1.setDimensions(0.12, 1);
		p1.add(new JLabel("Print Trees"));
		lO1.setDimensions(0.1, 1);
		p1.add(priBooBox);
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.737, 1);
		p2.add(new JLabel("Approximate Likelihood Ratio Test"));
		lO2.setDimensions(0.23, 1);
		p2.add(likRatBox);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == booStrBox) {
			if (booStrBox.getSelectedIndex() == 1) {
				likRatBox.setSelectedIndex(0);
				numRepField.setEnabled(true);
				numRepField.requestFocus();
				priBooBox.setEnabled(true);
			} else {
				numRepField.setEnabled(false);
				priBooBox.setEnabled(false);
			}
		} else if (e.getSource() == likRatBox) {
			if (likRatBox.getSelectedIndex() > 0) {
				numRepField.setEnabled(false);
				priBooBox.setEnabled(false);
				booStrBox.setSelectedIndex(0);
			}
		}
	}

	/**
	 * Class to customize a JTextField that the user is only able to type in
	 * integer values.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Class to customize a JTextField that the user is only able to type in
		 * integer values.
		 * 
		 * @param string
		 *            String : Default String as displayed when initialised.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isInteger(string)) {
				this.setText("");
			}
		}

		/**
		 * Tests whether the input value is convertible to an Integer value.
		 * 
		 * @param str
		 *            String : input String
		 * @return boolean : true if convertable, false otherwise.
		 */
		private boolean isInteger(String str) {
			try {
				Integer.parseInt(str);
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Class which extends PlainDocument for checking the input in the
		 * JTextField as the user types it in.
		 * 
		 * @author Christoph Knapp
		 * 
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isInteger(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}
	}

	/**
	 * Retrieves the number of replicates if user selected to use Bootstrapping
	 * and typed in a number in the text field.
	 * 
	 * @return String : number of replicates.
	 */
	public String getNumBootstrap() {
		if (booStrBox.getSelectedItem().toString().equals("yes")
				&& !numRepField.getText().equals("")) {
			return numRepField.getText();
		}
		return "";
	}

	public void focusGained(FocusEvent e) {
	}

	public void focusLost(FocusEvent e) {
		if (numRepField.getText().equals("")) {
			booStrBox.setSelectedIndex(0);
			priBooBox.setSelectedIndex(0);
			likRatBox.setSelectedIndex(1);
		}
	}
}
